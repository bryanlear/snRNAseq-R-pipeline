import csv
import gzip

MANIFEST = "config/samples.tsv"
HUMAN_H5 = "brain_data_human/GSE174367_snRNA-seq_filtered_feature_bc_matrix.h5"
HUMAN_METADATA = "brain_data_human/GSE174367_snRNA-seq_cell_meta.csv.gz"

#########################################################
##################### Human data ########################
#########################################################

with gzip.open(HUMAN_METADATA, "rt", newline="") as handle:
    HUMAN_SAMPLES = sorted({row["SampleID"] for row in csv.DictReader(handle)})

rule human_all:
    input:
        expand("results/GSE174367/02_qc/{sample}.rds", sample=HUMAN_SAMPLES),
        expand("results/GSE174367/02_qc/{sample}.thresholds.tsv", sample=HUMAN_SAMPLES)

rule import_human:
    input:
        h5=HUMAN_H5,
        metadata=HUMAN_METADATA,
        script="scripts/00_combined_h5.R"
    output:
        objects=expand("results/GSE174367/01_seurat/{sample}.rds", sample=HUMAN_SAMPLES),
        summary="results/GSE174367/01_seurat/import_summary.tsv",
        excluded="results/GSE174367/01_seurat/excluded_barcodes.tsv"
    params:
        output_dir="results/GSE174367/01_seurat",
        unmatched_policy=config.get("unmatched_policy", "exclude")
    threads: 1
    log:
        "logs/GSE174367/00_combined_h5.log"
    shell:
        """
        mkdir -p logs/GSE174367
        Rscript {input.script:q} {input.h5:q} {input.metadata:q} {params.output_dir:q} {params.unmatched_policy:q} > {log:q} 2>&1
        """

rule qc_human:
    input:
        rds="results/GSE174367/01_seurat/{sample}.rds",
        script="scripts/02_qc.R"
    output:
        rds="results/GSE174367/02_qc/{sample}.rds",
        summary="results/GSE174367/02_qc/{sample}.thresholds.tsv"
    params:
        nmads=3
    threads: 1
    log:
        "logs/GSE174367/02_qc/{sample}.log"
    shell:
        """
        mkdir -p logs/GSE174367/02_qc
        Rscript {input.script:q} {input.rds:q} {output.rds:q} {output.summary:q} {params.nmads} > {log:q} 2>&1
        """
        
#################### To call at any moment #################

rule human_mito_check_all:
    input:
        expand("results/GSE174367/02.1_mito_check/{sample}.cells.tsv", sample=HUMAN_SAMPLES),
        expand("results/GSE174367/02.1_mito_check/{sample}.summary.tsv", sample=HUMAN_SAMPLES)

rule qc_mito_summary:
    input:
        rds="results/GSE174367/02_qc/{sample}.rds",
        script="scripts/02.1_qc.R"
    output:
        cells="results/GSE174367/02.1_mito_check/{sample}.cells.tsv",
        summary="results/GSE174367/02.1_mito_check/{sample}.summary.tsv"
    params:
        mito_pattern="^MT-"
    threads: 1
    log:
        "logs/GSE174367/02.1_mito_check/{sample}.log"
    shell:
        """
        mkdir -p logs/GSE174367/02.1_mito_check
        Rscript {input.script:q} {input.rds:q} {output.cells:q} {output.summary:q} {params.mito_pattern:q} > {log:q} 2>&1
        """

##################################################################################################################
##################################################################################################################
##################################################################################################################
##################################################################################################################


#########################################################
############### Non-human primate data ##################
#########################################################

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