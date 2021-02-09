# RNA-seq analysis pipeline

The pipeline is build with Snakemake. It generates the alignments with STAR and quantifications with the updated Salmon algorithm.

## What do you need to run the pipeline?

* [`Snakemake`](https://snakemake.readthedocs.io/en/stable/) - tested on versions >= 5.11
* [`conda`](https://docs.conda.io/en/latest/) - I strongly recommend using conda, as it will automatically download and install all the necessary packages; it is possible to run the pipeline without it though. For the versions of the used software please refer to `utils/envs/` where you can browse environments specifications. You
* `fastq.gz` files (for now the pipeline accepts gzipped files) stored in the `reads/raw` directory
* `sample_sheet.tsv`: tab delimited file with sample ID (matching the fastq files names) in the first column
* `config.yaml`: there are some fields in the file that have to be filled accordingly to user specific set up.

This is how the directory tree will look like:

```bash
.
|-- Snakefile
|-- alignment
|   `-- quality
|-- logs
|   |-- MarkDuplicates
|   |-- STAR
|   |-- fastqc
|   |-- salmon
|   `-- trim_galore
|-- multiqc_data
|-- quants
|-- reads
|   |-- quality
|   |-- raw
|   `-- trimmed
|-- ref -> /path/to/ref
|-- sample_sheet.tsv -> sample_sheet_all.tsv
|-- sample_sheet_all.tsv
|-- sample_sheet_one.tsv
`-- utils
```

I tend to link `ref` to the reference path shared between projects not to copy the big files and refer to `ref` within the current directory. 

## Quick setup

Copy the content of this directory to the directory where the analysis will take place, make sure the reads are available in `reads/raw` subdirectory and, if you have `conda`, you can proceed.

The pipeline relies on the conda environments and it is strongly encouraged to use it. `Conda` to makes sure the pipeline is fully reproducible.

Example usage on the HPC.

```bash
snakemake --use-conda --profile cluster --jobs 50
```

It is recommended to first test the execution with `-n/--dry-run` to get the job counts and make sure there are no errors.

From experience, I would also recommend run the pipeline for one sample; make sure everything runs smoothly and then run the whole set. I do that by creating two sets of sample sheets and linking them to the expected `sample_sheet.tsv` accordingly. Example below:

```bash
head -n2 sample_sheet_all.tsv > sample_sheet_one.tsv
ln -s sample_sheet_one.tsv sample_sheet.tsv # make sure sample_sheet.tsv points to one sample only
snakemake -np --use-conda --profile cluster # dry run to make sure it works
snakemake --use-conda --profile cluster # run for one sample

unlink sample_sheet.tsv && ln -s sample_sheet_all.tsv sample_sheet.tsv
snakemake -np --use-conda --profile cluster # dry run to make sure it works for all samples
snakemake --use-conda --profile cluster
```

If something breaks - contact me! I'll do my best to help you out :)
