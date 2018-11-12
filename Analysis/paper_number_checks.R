##Abstract####
#Although 39% of murine protein-coding genes are essential for development
round(length(which(universe_df$lethal_mouse=="Y"))/length(which(universe_df$mouse_ko=="Y"))*100,0)

#we curate only 624 human disease genes (3% of protein-coding genes) linked to prenatal/infantile lethality
length(which(universe_df$human_lethal_B=="Y"))
round(length(which(universe_df$human_lethal_B=="Y"))/length(universe_df$gene)*100,0)

#then critically assessed features of 3,187 clinically relevant OMIM genes compared to 16,009 non-disease genes
length(which(universe_df$omim=="Y"))
length(which(is.na(universe_df$omim=="Y")))

#54% of all OMIM genes and 75% human lethal genes are linked to developmental lethality in knockout mice, compared to 34% non-disease genes
round(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))*100,0)
round(length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))*100,0)
round(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))*100,0)

#72.6% of OMIM genes are not constrained
round(length(which(universe_df$omim=="Y"&universe_df$constrained=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))*100,0)

#recessive disorders (only 8% constrained). 
round(length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="Y"))/
  length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))*100,0)

#We further curate 3423 ‘candidate developmental lethal’ human genes: essential for murine development or cellular viability, not yet linked to human disorders
length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|universe_df$cell_essential=="Y")))

##Introduction####
#Mouse Genome informatics (MGI) provides curated phenotype data for murine models with targeted knock-out of 6,991 protein-coding genes 
length(which(!is.na(universe_df$lethal_MGI)))

#IMPC no. genes
length(which(!is.na(universe_df$lethal_IMPC)))

# 75% known human lethal genes linked also to developmental lethal murine phenotypes
round(length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))*100,0)

##Results####
## Informatics datasets integrated into the Gene Discovery Toolkit ####
#Among 9397 genes for which murine phenotypic data was available, 39.2% result in pre-weaning lethality (< 3 weeks of age)
length(which(universe_df$mouse_ko=="Y"))
round(length(which(universe_df$lethal_mouse=="Y"))/length(which(universe_df$mouse_ko=="Y"))*100,1)

#assessing 15903/19196 (83%) of protein coding genes among eleven different cell lines 
length(which(universe_df$cell_ko=="Y"))
round(length(which(universe_df$cell_ko=="Y"))/length(universe_df$gene)*100,0)

#we define ‘2233 cell essential genes’ as those genes causing non-viability when knocked out in three or more cell lines 
length(which(universe_df$cell_essential=="Y"))

#Around 18% of human protein-coding genes are essential for cell viability in culture (aligning well with ~19% genes essential in yeast
round(length(which(universe_df$cell_essential_hits>=2))/length(which(universe_df$cell_ko=="Y"))*100,0)

#Only 21% of cell essential genes (478/2233) are linked currently to human disorders. 
length(which(universe_df$cell_essential=="Y"&universe_df$omim=="Y"))/length(which(universe_df$cell_essential=="Y"))

##Genetic constraint is a poor predictor of being a ‘disease gene’ - especially for recessive genes. ####
#Importantly, 72.6% of OMIM genes do not exhibit whole-gene genetic constraint
round(length(which(universe_df$omim=="Y"&universe_df$constrained=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))*100,1)

#54% of known clinically relevant OMIM genes are classified as tolerant to genetic variation 
round(length(which(universe_df$omim=="Y"&universe_df$any_constraint=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))*100,0)

#Levels of genetic constraint correlate with inheritance pattern, with human lethal genes showing highest levels of genetic constraint ####

#Autosomal recessive (AR) genes rarely demonstrate missense constraint (4.3%) or intolerance to LoF variation (8.4%). 
round(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))*100,1)
round(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))*100,1)

#genes associated with autosomal dominant (AD) disorders show significantly higher levels of missense constraint (29%), with 46% intolerant to LoF variation

round(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))*100,1)
round(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))*100,1)


#X-linked (XL) genes show intermediate levels of missense constraint (20.1% XL compared to 4.3% AR and 29% AD); though are the most intolerant of LoF variation (66.9%)
round(length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&universe_df$mis_z>=3.09))/length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&!is.na(universe_df$exac)))*100,1)
round(length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&universe_df$pLI>=0.9))/length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&!is.na(universe_df$exac)))*100,1)

#36 MT genes
length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))


#Manifestation of a severe-lethal animal phenotype correlates most strongly with being a disease gene####
#59% of OMIM genes linked to early lethality (developmental or neonatal) in knockout mice, compared to 37% non-disease genes
#among our curated list of 624 genes associated with human prenatal, perinatal or infantile lethality – 81% were also associated with a developmental (or neonatal) lethal murine phenotype
round(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))*100,0)
round(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))*100,0)
round(length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))*100,0)

data.frame(non_OMIM=c(length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y")),
                      length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y")),
                      length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))),
           OMIM=c(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y")),
                  length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y")),
                  length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))),
           Human_Lethal_B=c(length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y")),
                            length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y")),
                            length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="N"&!grepl("premature death",universe_df$all_MP_phen)))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))))

#our attention focusses intently on 2377 genes known essential for murine development, not yet linked to human disease
length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y")))

#Cell essential genes present strong candidates for unexplained infertility or early embryonic lethality####
#we noted 709 cell essential genes not known to be associated with human disorders, though linked to embryological lethality in mice 
length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="Y"))

#Unexpectedly, 54 genes classified as cell-essential were not linked to murine lethal phenotypes 
#X/54 genes were associated with sub-viable murine phenotypes and Y/54 with phenotypic abnormality of a certain cell type or organ 
length(which(universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"))
length(which(universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"&grepl("cellular phenotype",universe_df$high_MP_phen)))
length(which(universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"&grepl("premature death",universe_df$all_MP_phen)))



#Gene ontology analyses of the 3423 ‘candidate developmental lethal’ genes

ne<-read.table("Gene_lists/BINGO/candidates_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")

a<-ne$Genes.in.test.set[which(ne$Description%in%c("RNA binding","transcription regulator activity",
                                                  "transcription","transcription factor activity",
                                                  "DNA binding","nucleotide binding","translation",
                                                  "translation factor activity","nucleus"))]
dom<-unique(unlist(lapply(a, function(x) strsplit(as.character(x),"\\|"))))

length(which(universe_df$constrained[which(universe_df$gene%in%dom)]=="Y"))


b<-ne$Genes.in.test.set[which(ne$Description%in%c("mitochondrion"))]
rec<-unique(unlist(lapply(b, function(x) strsplit(as.character(x),"\\|"))))
length(which(universe_df$constrained[which(universe_df$gene%in%rec)]=="Y"))

