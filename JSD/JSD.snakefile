configfile: "/global/home/users/ryo10244201/sudmant_lab/gene_expression_variance/snakemake/config.json"

def get_all_inputs(wildcards):
    return expand(config['JSD'] + "{tissue}_JSD_result.csv", tissue=config["tissues_new"])

rule all:
    input:
        get_all_inputs

rule calc_JSD:
    input:
        config['path_new'] + "GTEx_Analysis-{tissue}.tsv"
    output:
        config['JSD'] + "{tissue}_JSD_result.csv"
    script:
        "../script/JSD/JSD_all.py"
