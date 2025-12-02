# RaqManzano/nf_rna_pipeline

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/RaqManzano/nf_rna_pipeline)
[![GitHub Actions CI Status](https://github.com/RaqManzano/nf_rna_pipeline/actions/workflows/nf-test.yml/badge.svg)](https://github.com/RaqManzano/nf_rna_pipeline/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/RaqManzano/nf_rna_pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/RaqManzano/nf_rna_pipeline/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/RaqManzano/nf_rna_pipeline)

## Introduction

**RaqManzano/nf_rna_pipeline** is a bioinformatics pipeline that process RNA-seq data. It will run:

1. QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) [default: false]
2. STAR alignment
3. salmon quantification
4. Genotyping with HaplotypeCaller
5. [`MultiQC`](http://multiqc.info/)

# Before you start

## Management of input/output files

Please also take a moment to think where the data will be downloaded.

As a suggestion keep your data tidy:

Create a folder e.g. “mkdir rna_analysis”, inside generate a folder with your sequencing id and the subfolders `data`, `metadata`, and `output`

```
rna_analysis
└──SequencingID1
   ├── data     # store your data (FQ/BAM)
   ├── metadata # store the info about the samples (e.g. samplesheet)
   └── output   # store the output for the analysis
```

Note: If using `clarity_tools` run clarity tools to download your data, please run it inside `data` folder this will generate a folder with the sequencing id that you can change name to `raw` with the following command

```
java -jar clarity-tools.jar -l <SLX-12345>
mv SLX1234 raw
```

This will download fastq files in a compressed format (“`fq.gz`”), quality control HTML reports which you can take a look at and csv file with sample information used for the sequencing (can be used as metadata).

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Generate samplesheet

The sample sheet should look like this:

```csv
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

or

```csv
sample,bam
CONTROL_REP1,AEG588A1_S1_L002_R1_001.bam
```

To generate a samplesheet from a list of fq.files specific for CRUK CI clarity-tools download:

```
python /path/to/pipeline/bin/make_input_file.py --path "/path/to/rna_analysis/SLX1234/data/raw/" --output SLX1234_samplesheet.txt
```

### Minimum params

The minimum set of parameters that should be provided are shown in the `yml` below.

```
input: '/path/to/samplesheet.csv'
outdir: '/path/to/outdir/'
read_length: 50 # please specify the read length of your reads for STAR alignment
transcriptome: '/path/to/transcriptome.fa'
gtf: '/path/to/reference.gtf'
fasta: '/path/to/genome.fa'
```

### Config

#### Env config

Please set up the config of your environment. For example, if you are running the pipeline in an HPC environment create a `hpc.config`:

```
singularity {
  enabled = true
  runOptions = "--no-home"
  pullTimeout = '60m' // Extend the timeout to 60 minutes
  autoMounts = true
}

process {
  executor          = 'slurm'
  clusterOptions    = "--account myaccount"
  queue             = "yourpartition"
  cache             = 'lenient'
  errorStrategy     = "retry"  // retry if error
  maxRetries = 2
}

executor { // these parameters avoid overcalling the scheduler
    queueSize         = 2005
    pollInterval      = '30 sec'
    queueStatInterval = '30 sec'
    submitRateLimit   = '50sec'
    exitReadTimeout   = "3 min"

}
```

### Usage in an HPC environment

It is highly recommended to run the pipeline in a terminal multiplexer such as `[screen](https://www.geeksforgeeks.org/linux-unix/screen-command-in-linux-with-examples/)` or `[tmux](https://github.com/tmux/tmux/wiki)`. A terminal multiplexer lets you open “virtual terminals” that continue running in the background, allowing processes to persist even if you close your window or disconnect from the HPC (think about it as a terminal window manager). This is important as the nextflow run by default needs to keep running to check jobs. By default `screen` is likely to be installed already in your environment.

To start a new session you can use the following commands:

```
# screen command:
screen -S demo
# tmux command:
tmux new -s demo
```

#### Tools config

A set of default config for STAR, salmon and HaplotypeCaller can be found in [config](conf/modules.config). For more information on the parameters of each tool please go to their docs.

## Basic command

```bash
nextflow run RaqManzano/nf_rna_pipeline \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

RaqManzano/nf_rna_pipeline was originally written by Raquel Manzano.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use RaqManzano/nf_rna_pipeline for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
