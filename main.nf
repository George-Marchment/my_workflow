#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * Tumour genomic profiling: from tumour sequencing reads to a clinical report.
 * Small variants are always called; the rest is switched on by the parameters.
 */

params.reads       = "data/*.fastq.gz"
params.reference   = "ref/GRCh38.fa"
params.hla_alleles = "HLA-A*02:01,HLA-B*07:02"
params.long_reads  = false   // Nanopore/PacBio instead of Illumina
params.structural  = true    // copy number and structural variants
params.neoantigens = true    // neoantigen prediction, for immunotherapy

// ---------- Short reads ----------

process FASTP {
    input:
    tuple val(id), path(reads)
    output:
    tuple val(id), path("${id}.trimmed.fastq.gz")
    script:
    """
    fastp -i ${reads} -o ${id}.trimmed.fastq.gz
    """
}

process BWA_MEM2 {
    input:
    tuple val(id), path(reads)
    output:
    tuple val(id), path("${id}.bam")
    script:
    """
    bwa-mem2 mem ${params.reference} ${reads} | samtools sort -o ${id}.bam
    """
}

// ---------- Long reads ----------

process NANOFILT {
    input:
    tuple val(id), path(reads)
    output:
    tuple val(id), path("${id}.filtered.fastq.gz")
    script:
    """
    gunzip -c ${reads} | NanoFilt -q 10 | gzip > ${id}.filtered.fastq.gz
    """
}

process MINIMAP2 {
    input:
    tuple val(id), path(reads)
    output:
    tuple val(id), path("${id}.bam")
    script:
    """
    minimap2 -ax map-ont ${params.reference} ${reads} | samtools sort -o ${id}.bam
    """
}

// ---------- Preprocessing ----------

process MARK_DUPLICATES {
    input:
    tuple val(id), path(bam)
    output:
    tuple val(id), path("${id}.dedup.bam")
    script:
    """
    gatk MarkDuplicates -I ${bam} -O ${id}.dedup.bam -M ${id}.metrics.txt
    """
}

// ---------- Small variants ----------

process MUTECT2 {
    input:
    tuple val(id), path(bam)
    output:
    tuple val(id), path("${id}.mutect2.vcf.gz"), path("${id}.f1r2.tar.gz")
    script:
    """
    gatk Mutect2 -R ${params.reference} -I ${bam} --f1r2-tar-gz ${id}.f1r2.tar.gz -O ${id}.mutect2.vcf.gz
    """
}

process LEARN_ORIENTATION {
    input:
    tuple val(id), path(vcf), path(f1r2)
    output:
    tuple val(id), path("${id}.orientation.tar.gz")
    script:
    """
    gatk LearnReadOrientationModel -I ${f1r2} -O ${id}.orientation.tar.gz
    """
}

process FILTER_CALLS {
    input:
    tuple val(id), path(vcf), path(f1r2), path(orientation)
    output:
    tuple val(id), path("${id}.filtered.vcf.gz")
    script:
    """
    gatk FilterMutectCalls -R ${params.reference} -V ${vcf} --ob-priors ${orientation} -O ${id}.filtered.vcf.gz
    """
}

process STRELKA2 {
    input:
    tuple val(id), path(bam)
    output:
    tuple val(id), path("${id}.strelka2.vcf.gz")
    script:
    """
    configureStrelkaSomaticWorkflow.py --tumorBam ${bam} --referenceFasta ${params.reference} --runDir strelka
    strelka/runWorkflow.py -m local
    mv strelka/results/variants/somatic.snvs.vcf.gz ${id}.strelka2.vcf.gz
    """
}

process VEP {
    input:
    tuple val(id), path(mutect2), path(strelka2)
    output:
    tuple val(id), path("${id}.annotated.vcf")
    script:
    """
    bcftools concat -a ${mutect2} ${strelka2} | vep --vcf --cache --output_file ${id}.annotated.vcf
    """
}

// ---------- Neoantigens ----------

process NETMHCPAN {
    input:
    tuple val(id), path(vcf)
    output:
    tuple val(id), path("${id}.netmhcpan.tsv")
    script:
    """
    pvacseq run ${vcf} ${id} ${params.hla_alleles} NetMHCpan . -e1 9
    mv MHC_Class_I/${id}.all_epitopes.tsv ${id}.netmhcpan.tsv
    """
}

process MHCFLURRY {
    input:
    tuple val(id), path(vcf)
    output:
    tuple val(id), path("${id}.mhcflurry.tsv")
    script:
    """
    pvacseq run ${vcf} ${id} ${params.hla_alleles} MHCflurry . -e1 9
    mv MHC_Class_I/${id}.all_epitopes.tsv ${id}.mhcflurry.tsv
    """
}

process RANK_NEOANTIGENS {
    input:
    tuple val(id), path(netmhcpan), path(mhcflurry)
    output:
    tuple val(id), path("${id}.neoantigens.tsv")
    script:
    """
    rank_neoantigens.py ${netmhcpan} ${mhcflurry} > ${id}.neoantigens.tsv
    """
}

// ---------- Structural variants ----------

process CNVKIT {
    input:
    tuple val(id), path(bam)
    output:
    tuple val(id), path("${id}.cns")
    script:
    """
    cnvkit.py batch ${bam} --fasta ${params.reference} --output-dir .
    """
}

process MANTA {
    input:
    tuple val(id), path(bam)
    output:
    tuple val(id), path("${id}.sv.vcf.gz")
    script:
    """
    configManta.py --tumorBam ${bam} --referenceFasta ${params.reference} --runDir manta
    manta/runWorkflow.py
    mv manta/results/variants/tumorSV.vcf.gz ${id}.sv.vcf.gz
    """
}

process PURPLE {
    input:
    tuple val(id), path(cns), path(sv)
    output:
    tuple val(id), path("${id}.purity.tsv")
    script:
    """
    purple -tumor ${id} -somatic_sv_vcf ${sv} -copy_number ${cns} -output_dir .
    """
}

// ---------- Report ----------

process TUMOUR_REPORT {
    input:
    path(results)
    output:
    path("report.html")
    script:
    """
    build_report.py ${results} > report.html
    """
}

// ---------- Subworkflows ----------

workflow SHORT_READS {
    take:
    reads
    main:
    FASTP(reads)
    BWA_MEM2(FASTP.out)
    emit:
    bam = BWA_MEM2.out
}

workflow LONG_READS {
    take:
    reads
    main:
    NANOFILT(reads)
    MINIMAP2(NANOFILT.out)
    emit:
    bam = MINIMAP2.out
}

workflow PREPROCESSING {
    take:
    reads
    main:
    if (params.long_reads) {
        LONG_READS(reads)
        bam = LONG_READS.out.bam
    } else {
        SHORT_READS(reads)
        bam = SHORT_READS.out.bam
    }
    MARK_DUPLICATES(bam)
    emit:
    bam = MARK_DUPLICATES.out
}

workflow VARIANT_CALLING {
    take:
    bam
    main:
    MUTECT2(bam)
    LEARN_ORIENTATION(MUTECT2.out)
    FILTER_CALLS(MUTECT2.out.join(LEARN_ORIENTATION.out))
    STRELKA2(bam)
    emit:
    vcfs = FILTER_CALLS.out.join(STRELKA2.out)
}

workflow NEOANTIGEN_PREDICTION {
    take:
    vcf
    main:
    NETMHCPAN(vcf)
    MHCFLURRY(vcf)
    RANK_NEOANTIGENS(NETMHCPAN.out.join(MHCFLURRY.out))
    emit:
    ranked = RANK_NEOANTIGENS.out
}

workflow SMALL_VARIANTS {
    take:
    bam
    main:
    VARIANT_CALLING(bam)
    VEP(VARIANT_CALLING.out.vcfs)
    if (params.neoantigens) {
        NEOANTIGEN_PREDICTION(VEP.out)
        results = VEP.out.mix(NEOANTIGEN_PREDICTION.out.ranked)
    } else {
        results = VEP.out
    }
    emit:
    results
}

workflow STRUCTURAL_VARIANTS {
    take:
    bam
    main:
    CNVKIT(bam)
    MANTA(bam)
    PURPLE(CNVKIT.out.join(MANTA.out))
    emit:
    purity = PURPLE.out
}

workflow GENOMIC_PROFILING {
    take:
    bam
    main:
    SMALL_VARIANTS(bam)
    if (params.structural) {
        STRUCTURAL_VARIANTS(bam)
        results = SMALL_VARIANTS.out.results.mix(STRUCTURAL_VARIANTS.out.purity)
    } else {
        results = SMALL_VARIANTS.out.results
    }
    emit:
    results
}

workflow {
    reads = Channel.fromPath(params.reads).map { file -> tuple(file.simpleName, file) }
    PREPROCESSING(reads)
    GENOMIC_PROFILING(PREPROCESSING.out.bam)
    TUMOUR_REPORT(GENOMIC_PROFILING.out.results.map { id, file -> file }.collect())
}
