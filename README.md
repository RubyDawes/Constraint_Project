# Constraint_Project
Integration of available data from OMIM (mendelian disease), ExAC (genetic constraint), MGI (lethal mouse knockouts), CRISPR screens (cell essentiality) for analysis

## Data Sources
### 1. OMIM
mim2gene.txt downloaded from [OMIM Data Downloads page](https://www.omim.org/downloads/) 
genemap2.txt downloaded with license on 2018-03-10
- file location: Gene_lists/OMIM/...


### 2. ExAC
from [Lek et al., Analysis of protein-coding genetic variation in 60,706 humans, 2016](https://www.nature.com/articles/nature19057)

ExAC scores downloaded from [ExAC Downloads page](http://exac.broadinstitute.org/downloads) 
- file location: Gene_lists/ExAC/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt


### 3. Mouse knockouts
HMD_HumanPhenotype, MGI_PhenoGenoMP, MGI_PhenotypicAllele downloaded from [MGI Data and Statistical Reports page](http://www.informatics.jax.org/downloads/reports/index.html) and converted to xlsx for ease of reading in R/excel
- filt location: Gene_lists/MGI_MPs/...

Gene_lists/MGI_MPs/MGI_lethal_phenotypes is a result of manual curation of MP terms to find lethal MP terms

### 4. CRISPR screens
- Gene_lists/Cell_KOs/blomen_essentials
from [Blomen et al., Gene essentiality and synthetic lethality in haploid human cells, 2015](http://science.sciencemag.org/content/350/6264/1092) tables S1 and S2. essentialome from table s3

- Gene_lists/Cell_KOs/hart_essentials
from [Hart et al., High-Resolution CRISPR Screens Reveal Fitness Genes and Genotype-Specific Cancer Liabilities, 2015](https://www.cell.com/cell/fulltext/S0092-8674(15)01495-6) table S2 [bayes cutoffs found in supplementary experimental procedures]

- Gene_lists/Cell_KOs/wang_essentials
from [Wang et al., Identification and characterization of essential genes in the human genome, 2015](http://science.sciencemag.org/content/350/6264/1096) table S3 and S4


### 5. Universe
list of all protein-coding genes downloaded from [HGNC statistics & Downloads page](https://www.genenames.org/cgi-bin/statistics)
- file location: Gene_lists/Universe/gene_with_protein_product.txt

## Running the scripts
'Master.R' describes each script and guides you through the order in which to run scripts, and whether to load from rda to avoid slow scripts


