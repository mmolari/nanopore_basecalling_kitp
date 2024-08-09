import pathlib
import numpy as np

# create log directory
log_fld = pathlib.Path("log")
log_fld.mkdir(exist_ok=True)

kit = config["kit"]
no_barcoding_kits = ["SQK-LSK114"]

# pod5 files numbers:
pod5_fld = "nanopore_runs/" + config["run_id"] + "/pod5"
pod5_ids = glob_wildcards(pod5_fld + "/{sample}.pod5").sample
batch_size = 10
N_batches = int(np.ceil(len(pod5_ids) / batch_size))

print(f"pod5 folder: {pod5_fld}")
print(f"samples: {pod5_ids}")
print(f"batch size: {batch_size}")
print(f"N batches: {N_batches}")


rule download_dorado_model:
    output:
        directory("dorado_models/{model}"),
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} download --model {wildcards.model} --directory dorado_models
        """


def bath_n_to_batch_files(wildcards):
    n = int(wildcards.n)
    pod5_selected = pod5_ids[n * batch_size : (n + 1) * batch_size]
    return expand(pod5_fld + "/{sample}.pod5", sample=pod5_selected)


rule create_batch:
    input:
        bath_n_to_batch_files,
    output:
        batch=directory("nanopore_runs/{run_id}/batches/batch_{n}"),
    shell:
        """
        mkdir -p {output}
        # create symlinks
        for f in {input}; do
            ln -s ../../pod5/$(basename $f) {output}
        done
        """


rule basecall:
    input:
        rds=rules.create_batch.output.batch,
        mdl=expand(rules.download_dorado_model.output, model=config["dorado_model"]),
    output:
        bam="basecalled/{run_id}/basecalled_bam/batch_{n}.bam",
    params:
        kit=lambda w: "" if kit in no_barcoding_kits else f"--kit-name {kit}",
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} basecaller {input.mdl} {input.rds} {params.kit} > {output}
        """


checkpoint demux:
    input:
        bam=rules.basecall.output,
    output:
        bcd=directory("basecalled/{run_id}/barcodes_bam/batch_{n}"),
    params:
        dorado_bin=config["dorado_bin"],
    shell:
        """
        {params.dorado_bin} demux {input.bam} --no-classify --output-dir {output}
        """


def all_barcodes_bam(wildcards):
    files = []
    for n_batch in range(N_batches):
        # is barcode in batch output?
        batch_bam_fld = pathlib.Path(
            checkpoints.demux.get(run_id=config["run_id"], n=n_batch).output["bcd"]
        )
        batch_bam_file = batch_bam_fld / f"{wildcards.barcode}.bam"
        print("batch bam fld:", batch_bam_fld)
        print("batch bam file:", batch_bam_file)
        print("exists:", batch_bam_file.exists())

        if (batch_bam_fld / f"{wildcards.barcode}.bam").exists():
            files += expand(
                "basecalled/{run_id}/barcodes_bam/batch_{n}/{barcode}.bam",
                barcode=wildcards.barcode,
                run_id=config["run_id"],
                n=n_batch,
            )
    print(f"requested files for {wildcards.barcode}:\n", files)
    return files


rule collect_fastq:
    input:
        all_barcodes_bam,
    output:
        "basecalled/{run_id}/reads/{barcode}.fastq.gz",
    shell:
        """
        for f in {input}; do
            samtools fastq $f | gzip >> {output}
        done
        """


rule to_fastq_no_barcoding:
    input:
        rules.basecall.output.bam,
    output:
        "basecalled/{run_id}/basecalled_bam/batch_{n}.fastq.gz",
    shell:
        """
        samtools fastq {input} | gzip > {output}
        """


rule collect_fastq_no_barcoding:
    input:
        expand(
            rules.to_fastq_no_barcoding.output,
            run_id=config["run_id"],
            n=range(N_batches),
        ),
    output:
        "basecalled/{run_id}/reads.fastq.gz",
    shell:
        """
        cat {input} > {output}
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


# rule summary:
#     input:
#         bam=rules.basecall.output.bam,
#     output:
#         "basecalled/{run_id}/summary.tsv",
#     params:
#         dorado_bin=config["dorado_bin"],
#     shell:
#         """
#         {params.dorado_bin} summary {input.bam} > {output}
#         """


def all_fastq(wildcards):
    if kit in no_barcoding_kits:
        # return without demultiplexing
        return expand(rules.collect_fastq_no_barcoding.output, run_id=config["run_id"])
    else:
        all_barcodes = []
        for n_batch in range(N_batches):
            B = pathlib.Path(
                checkpoints.demux.get(run_id=config["run_id"], n=n_batch).output["bcd"]
            ).glob("*.bam")
            all_barcodes.extend([x.stem for x in B])
        all_barcodes = list(set(all_barcodes))
        # return all fastq files
        return expand(
            rules.collect_fastq.output, barcode=all_barcodes, run_id=config["run_id"]
        )


rule all:
    input:
        # expand(rules.summary.output, run_id=config["run_id"]),
        expand(rules.config_info.output, run_id=config["run_id"]),
        all_fastq,


localrules:
    download_dorado_model,
    create_batch,
