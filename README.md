# Gene Expression Aging V1.0
This repository contains all code for "Tissue-specific impacts of aging and genetics on gene expression patterns in humans" by Yamamoto and Chung et al. 2022

Abstract
--------


Age is the primary risk factor for many common human diseases including heart disease, Alzheimer’s dementias, cancers, and diabetes. Determining how and why tissues age differently is key to understanding the onset and progression of such pathologies. Here, we set out to quantify the relative contributions of genetics and aging to gene expression patterns from data collected across 27 tissues from 948 humans. We show that age impacts the predictive power of expression quantitative trait loci across several tissues. Jointly modelling the contributions of age and genetics to transcript level variation we find that the heritability (<img src="https://render.githubusercontent.com/render/math?math=h^2">) of gene expression is largely consistent among tissues. In contrast, the average contribution of aging to gene expression variance varied by more than 20-fold among tissues with (<img src="https://render.githubusercontent.com/render/math?math=R^2_{age} > h^2">) in 5 tissues. We find that the coordinated decline of mitochondrial and translation factors is a widespread signature of aging across tissues. Finally, we show that while in general the force of purifying selection is stronger on genes expressed early in life compared to late in life as predicted by Medawar’s hypothesis, a handful of highly proliferative tissues exhibit the opposite pattern. These non-Medawarian tissues exhibit high rates of cancer and age-of-expression associated somatic mutations in cancer. In contrast, gene expression variation that is under genetic control is strongly enriched for genes under relaxed constraint. Together we present a novel framework for predicting gene expression phenotypes from genetics and age and provide insights into the tissue-specific relative contributions of genes and the environment to phenotypes of aging.

Folder Structure
---------------
- **eQTL:** code used to generate the datasets in eQTL analysis
- **JSD:** code used to generate the datasets in JSD analysis
- **Multi-snp:** code used to generate the datasets in Multi-snp model analysis
- **PEER:** folder contains PEER factors corrected for age


Installation
---------------
- Download GTEx v8 from https://www.gtexportal.org/home/datasets
- Download complete results of model from https://doi.org/10.5281/zenodo.6533954

Authors
-------

* **Ryo Yamamoto** <a itemprop="sameAs" content="https://orcid.org/0000-0003-3134-145X" href="https://orcid.org/0000-0003-3134-145X" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> - *(Co-Lead Author) All analyses, manuscript* - [ryo1024](https://github.com/ryo1024)
* **Ryan Chung** - *(Co-Lead Author) All analyses, manuscript*
* **Juan M Vazquez** <a itemprop="sameAs" content="https://orcid.org/0000-0001-8341-2390" href="https://orcid.org/0000-0001-8341-2390" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> - *Manuscript* - [docmanny](https://vazquez.bio)
* **Huanjie Sheng** - *Manuscript*
* **Nilah Ioannidis** - *(Coresponding Author) Manuscript*
* **Peter Sudmant** <a itemprop="sameAs" content="https://orcid.org/0000-0002-9573-8248" href="https://orcid.org/0000-0002-9573-8248" target="orcid.widget" rel="me noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> - *(Coresponding Author) Manuscript* - [petersudmant](https://github.com/petersudmant) [Lab Website](https://www.sudmantlab.org)

