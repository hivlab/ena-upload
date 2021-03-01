import pandas as pd
from snakemake.utils import validate
import os
import re

configfile: "data/config.yaml"

WRAPPER_PREFIX="https://raw.githubusercontent.com/avilab/virome-wrappers"
HOST_GENOME="/Users/taavi/Projects/databases/human_g1k_v37.fasta"

with open(config["batches"], "r") as h:
    batches = h.read().splitlines()
RUNS=list(set([re.sub("_R[1,2]_001.fastq.gz", "", e.split("/")[-1:][0]) for e in batches]))


def get_reads(wildcards):
    pair=[e for e in batches if wildcards.run in e]
    input={}
    for e in pair:
        input.update({"in1": e} if re.search("_R1", e) else {"in2": e})
    return input


rule all:
    input:
        "output/data.txt",
        "output/studies.tsv",
        "output/samples.tsv",
        "output/experiments.tsv",
        "output/runs.tsv",
        expand("output/data/{run}_R{pair}_clean.fastq.gz", run = RUNS, pair = [1, 2])


rule form:
    input:
        samples="data/sample_metadata.csv",
        batches="data/batches.txt",
        config="data/config.yaml",
    output:
        xlsx="output/ena_submission.xlsx",
        data="output/data.txt"
    log: 
        "logs/form.log"
    shell:
        """
        (Rscript scripts/submission-metadata-excel.R -s {input.samples} -b {input.batches} -c {input.config} -x {output.xlsx} -d {output.data}) 2> {log}
        """


rule parse_form:
    input:
        rules.form.output.xlsx
    output:
        "output/studies.tsv",
        "output/samples.tsv",
        "output/experiments.tsv",
        "output/runs.tsv",
    log: 
        "logs/parse_form.log"
    params:
        extra="--action add --vir"
    shell:
        """
        (python scripts/process_xlsx.py --form {input[0]} --out_dir output {params.extra}) 2> {log}
        """


# Quality trimming
rule trim:
    input:
        unpack(get_reads)
    output:
        out1=temp("output/data/{run}_R1_trimmed.fastq.gz"),
        out2=temp("output/data/{run}_R2_trimmed.fastq.gz"),
    log:
        "logs/{run}_trim.log",
    params:
        extra=(
            "ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 ordered"
        ),
    resources:
        runtime=lambda wildcards, attempt: 180 + (attempt * 60),
        mem_mb=16000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbduk"



# Remove host sequences
rule maphost:
    input:
        in1=rules.trim.output.out1,
        in2=rules.trim.output.out2,
        ref=HOST_GENOME,
    output:
        outu=temp("output/data/{run}_interleaved_clean.fastq.gz"),
        statsfile="output/data/{run}_maphost.txt",
    log:
        "logs/{run}_maphost.log",
    params:
        extra="minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 qtrim=rl trimq=10 untrim usemodulo",
    resources:
        runtime=lambda wildcards, attempt: attempt * 360,
        mem_mb=16000,
    threads: 8
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/bbwrap"


rule reformat:
    input:
        rules.maphost.output.outu,
    output:
        out="output/data/{run}_R1_clean.fastq.gz",
        out2="output/data/{run}_R2_clean.fastq.gz", 
    log:
        "logs/{run}_reformat.log",
    params:
        extra="",
    resources:
        runtime=360,
        mem_mb=4000,
    wrapper:
        f"{WRAPPER_PREFIX}/v0.6/bbtools/reformat"
