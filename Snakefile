import pathlib

# create log directory
log_fld = pathlib.Path("log")
log_fld.mkdir(exist_ok=True)


rule download_dorado_model:
    output:
        directory("dorado_models/{model}"),
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} download --model {wildcards.model} --directory dorado_models
        """


rule basecall:
    input:
        rds="nanopore_runs/{run_id}/pod5",
        mdl=expand(rules.download_dorado_model.output, model=config["dorado_model"]),
    output:
        bam="basecalled/{run_id}/basecalled.bam",
    params:
        kit=config["kit"],
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} basecaller {input.mdl} {input.rds} --kit-name {params.kit} > {output}
        """


checkpoint demux:
    input:
        bam=rules.basecall.output,
    output:
        bcd=directory("basecalled/{run_id}/barcodes_bam"),
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} demux {input.bam} --no-classify --output-dir {output}
        """


rule to_fastq:
    input:
        "basecalled/{run_id}/barcodes_bam/{barcode}.bam",
    output:
        "basecalled/{run_id}/barcodes_fastq/{barcode}.fastq.gz",
    shell:
        """
        samtools fastq {input} | gzip > {output}
        """


rule config_info:
    output:
        "basecalled/{run_id}/run_info.txt",
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        echo "{config}" > {output}
        sed -i 's/, /\\n/g' {output}
        # log git commit
        echo 'git commit: ' >> {output}
        git rev-parse HEAD >> {output}
        # log dorado version
        echo 'dorado version: ' >> {output}
        {params.dorado_bin} --version 2>> {output}
        # log samtools version
        echo 'samtools version: ' >> {output}
        samtools --version | head -n1 >> {output}
        """


rule summary:
    input:
        bam=rules.basecall.output.bam,
    output:
        "basecalled/{run_id}/summary.tsv",
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} summary {input.bam} > {output}
        """


def all_fastq(wildcards):
    all_barcodes = pathlib.Path(
        checkpoints.demux.get(run_id=config["run_id"]).output["bcd"]
    ).glob("*.bam")
    # strip the suffix
    all_barcodes = [x.stem for x in all_barcodes]
    return expand(rules.to_fastq.output, barcode=all_barcodes, run_id=config["run_id"])


rule all:
    input:
        expand(rules.demux.output, run_id=config["run_id"]),
        expand(rules.summary.output, run_id=config["run_id"]),
        expand(rules.config_info.output, run_id=config["run_id"]),
        all_fastq,


localrules:
    download_dorado_model,
