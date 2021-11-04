configfile: "/global/home/users/ryo10244201/sudmant_lab/gene_expression_variance/snakemake/config.json"

def get_all_inputs(wildcards):
    return expand(config['JSD_eQTL_new'] + "{tissue}_lm_results_new.tsv", tissue=config["tissues_gt50"])

rule all:
    input:
        get_all_inputs

rule calc_JSD:
    input:
        config['eQTL'] + "GTEx_Analysis_v8_eQTL_expression_matrices/{tissue}.v8.normalized_expression.bed.gz"
    output:
        config['JSD_eQTL_new'] + "{tissue}_lm_results_new.tsv"
    script:
        "../script/JSD/JSD_eQTL.py"
