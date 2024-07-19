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
        rds=directory("nanopore_runs/{run_id}/pod5"),
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


rule demux:
    input:
        bam=rules.basecall.output,
    output:
        directory("basecalled/{run_id}/barcodes"),
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} demux {input.bam} --no-classify --output-dir {output}
        """


rule all:
    input:
        expand(rules.demux.output, run_id=config["run_id"]),


rule clear:
    shell:
        """
        rm log/*
        """


localrules:
    download_dorado_model,
