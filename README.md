# nanopore basecalling with dorado

This repositoray contains code for the nanopore basecalling of data collected during 2024 qbio course at KITP.

> [!WARNING]  
> This repository is a work in progress.

## setup

- download [dorado](https://github.com/nanoporetech/dorado) binaries, and in the snakefile set `DORADO` equal to the binary path.
- create an environment with [snakemake](https://snakemake.readthedocs.io/en/stable/) and the [slurm plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html).