# nanopore basecalling with dorado

This repositoray contains code for the nanopore basecalling of data collected during 2024 qbio course at KITP.

> [!WARNING]  
> This repository is a work in progress.

## setup

- download [dorado](https://github.com/nanoporetech/dorado) binaries, and in the snakefile set `DORADO` equal to the binary path.
- create an environment with [snakemake](https://snakemake.readthedocs.io/en/stable/) v.7 and samtools.
    ```sh
    conda create -n snakemake -c conda-forge -c bioconda snakemake=7 samtools
    ```

## data preparation

- upload the `.pod5` files in a folder `nanopore_runs/<run-id>/pod5`
- prepare in the folder a `run_config.yml` file, using [this as template](nanopore_runs/example_run/run_config.yml)

## running the pipeline

after activating the conda environment, run the sequencing with the command

```sh
snakemake all \
    --profile profiles/cluster \
    --configfile nanopore_runs/<run-id>/run_config.yml
```

## output

basecalled reads will be saved in the `basecalled/<run-id>/` folder.