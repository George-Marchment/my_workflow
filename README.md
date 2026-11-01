# Tumour genomic profiling (demo)

A demo [Nextflow](https://www.nextflow.io/) workflow that takes tumour sequencing reads and produces a single HTML report describing each tumour's mutations.

It is meant to show what such a pipeline looks like, not to be used on real data. See [Limitations](#limitations).

## What it does



| Stage | What it answers | Tools |
|---|---|---|
| Preprocessing, short reads | Where does each read belong in the genome? | fastp, BWA-MEM2, GATK MarkDuplicates |
| Preprocessing, long reads | Same, for Nanopore data | NanoFilt, minimap2, GATK MarkDuplicates |
| Small variants | Which small mutations does the tumour carry, and what do they affect? | Mutect2, Strelka2, VEP |
| Neoantigens | Which mutations could produce a target for immunotherapy? | pVACseq with NetMHCpan and MHCflurry, then a ranking step |
| Structural variants | Which large regions are gained, lost or rearranged, and how pure is the sample? | CNVkit, Manta, PURPLE |
| Report | What are the combined results? | `build_report.py` |

Small variants are always called. The other analyses are switched on or off with the parameters below.

## Parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--reads` | `data/*.fastq.gz` | Input reads, one FASTQ file per sample. The sample name comes from the file name. |
| `--reference` | `ref/GRCh38.fa` | Reference genome |
| `--hla_alleles` | `HLA-A*02:01,HLA-B*07:02` | HLA alleles used for neoantigen prediction |
| `--long_reads` | `false` | Use the long-read branch instead of the short-read one |
| `--structural` | `true` | Run the structural variant analysis |
| `--neoantigens` | `true` | Run the neoantigen prediction |

## Usage

```bash
# Everything, on Illumina reads
nextflow run main.nf --reads 'data/*.fastq.gz' --reference ref/GRCh38.fa

# Long reads, small variants only
nextflow run main.nf --long_reads true --structural false --neoantigens false
```

## Output

The workflow produces one `report.html` covering all samples. It has no `publishDir`, so the report stays in Nextflow's `work/` directory.

## Limitations

- **Helper scripts are missing.** `rank_neoantigens.py` and `build_report.py` are called but not included, so a full run fails without them.
- **Tools are not packaged.** There is no container or Conda environment, so every tool must already be on the `PATH`.
- **Tumour only.** There is no matched normal sample, so inherited variants are not separated from tumour mutations as they would be in a clinical pipeline.
- **One FASTQ per sample.** Paired-end reads are not handled.
- **One HLA type for all samples.** The same alleles are applied to every sample rather than being typed per patient.
- **Nanopore settings only.** The long-read branch uses minimap2's `map-ont` preset, so PacBio reads would need a different one.
