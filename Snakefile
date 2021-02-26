
rule all:
    input:
        "output/data.txt",
        "output/studies.tsv",
        "output/samples.tsv",
        "output/experiments.tsv",
        "output/runs.tsv",


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
