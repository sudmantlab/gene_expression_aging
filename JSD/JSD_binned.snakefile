configfile: "/global/home/users/ryo10244201/sudmant_lab/variance_tissue/snakemake/config.json"

def get_all_inputs(wildcards):
    return expand(config['JSD_binned'] + "{tissue}_dist.csv", tissue=config["tissues"])

rule all:
    input:
        get_all_inputs

rule calc_JSD:
    input:
        config['path_new'] + "GTEx_Analysis-{tissue}.tsv"
    output:
        config['JSD_binned'] + "{tissue}_dist.csv"
    script:
        "../script/JSD_binned.py"
