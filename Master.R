#clears variables, loads dependencies, loads updated hgnc symbol table (for checking gene symbols), initialises universe_df 
##which will be added to with info from other scripts, loads universe which will be used as reference list of 'all protein-coding genes'
source("setup.R")

#functions for making pretty graphs
source("plot_functions.R")

#get omim df with filtered list of OMIM genes with mim numbers, phenotypes, and inheritance
##criteria for high confidence disease genes: phenotype mapping key (3) "molecular basis is known", filter out any nondiseases, 
###susceptibilities, provisional links, somatic mutations, removing non protein-coding genes (noncoding RNAs, complex loci etc)
source("1.OMIM_filtering/OMIM_parse_and_filter.R")

#get list of NMD genes
source("1.OMIM_filtering/nmd_genes.R")

#get exac scores for protein coding genes
source("2.ExAC_constraint/exac_constraint.R")

#take all genes with mouse knockout/hypomorph phenotypes from MGI, see which ones are lethal and what the lethal phenotypes are
##for list of lethal phenotypes see Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx
###warning: pretty time consuming (takes about ~1.5 mins to run) recommend loading from output spreadsheet or rda
#source("3.Mouse_phenotypes/MGI_parse.R")
#mgi <- read.xlsx("output/spreadsheets/MGI_genes_with_phenotypes.xlsx)
load("output/Data/mgi.rda")

#df of genes tested in wang, blomen and hart CRISPR knockout papers and whether they were lethal or not in each of the 11 cell lines
##tested throughout the three papers
#source("4.Cell_Knockouts/cell_essential_genes.R")
load("output/Data/cell_KOs.rda")

#puts all info on disease status, phenotype, exac scores, mouse phenotypes and cell knockouts for each protein-coding gene in HGNC into one 
##dataframe 'universe'
source("Universe/making_universe_df.R")
load("output/Data/universe_df.rda")

