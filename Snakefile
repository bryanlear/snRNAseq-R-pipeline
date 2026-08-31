import csv
MANIFEST = "config/samples.tsv"

with open(MANIFEST, newline="") as handle:
    sample_rows = list(csv.DictReader(handle, delimiter="\t"))

SAMPLES = {
    row["sample_id"]: row
    for row in sample_rows}
SAMPLE_IDS = list(SAMPLES)

rule all:
    input:
        expand("results/01_seurat/{sample}.rds",sample=SAMPLE_IDS)

rule create_seurat:
    input:
        h5=lambda wildcards: SAMPLES[wildcards.sample]["h5_file"],
        script="scripts/01_create_seurat.R"
    output:
        rds="results/01_seurat/{sample}.rds"
    params:
        condition=lambda wildcards: SAMPLES[wildcards.sample]["condition"]
    threads: 1
    log:
        "logs/01_create_seurat/{sample}.log"
    shell:
        """
        mkdir -p logs/01_create_seurat

        Rscript {input.script:q} \
            {wildcards.sample:q} \
            {params.condition:q} \
            {input.h5:q} \
            {output.rds:q} \
            > {log:q} 2>&1
        """