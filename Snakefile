DORADO = "~/apps/dorado-0.7.2-linux-x64/bin/dorado"


rule download_dorado_model:
    output:
        directory("dorado_models/{model}"),
    params:
        dorado_bin=DORADO,
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
        dorado_bin=DORADO,
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
        dorado_bin=DORADO,
    shell:
        """
        {params.dorado_bin} demux {input.bam} --no-classify --output-dir {output}
        """


rule all:
    input:
        expand(rules.demux.output, run_id=config["run_id"]),


localrules:
    download_dorado_model,
