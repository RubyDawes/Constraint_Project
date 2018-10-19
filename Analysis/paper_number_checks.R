#54% of all OMIM genes and 75% human lethal genes are linked to developmental lethality in knockout mice, compared to 34% non-disease genes
length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))
length(which(universe_df$human_lethal_B=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$human_lethal_B=="Y"&universe_df$mouse_ko=="Y"))
length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))

#72.6% of OMIM genes are not constrained
length(which(universe_df$omim=="Y"&universe_df$constrained=="N"))/length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))

#recessive disorders (only 8% constrained). 
length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="Y"))/
  length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))

#We further curate 3423 ‘candidate developmental lethal’ human genes: essential for murine development or cellular viability, not yet linked to human disorders
length(which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|universe_df$cell_essential=="Y")))


#we identified only 624 genes (3% of human protein coding genes) linked currently to prenatal or infantile lethality identified
length(which(universe_df$human_lethal_B=="Y"))
length(which(universe_df$human_lethal_B=="Y"))/length(universe_df$gene)*100

#39 % of murine genes are essential for murine development
length(which(universe_df$lethal_mouse=="Y"))/length(which(universe_df$mouse_ko=="Y"))

#Among 9397 genes for which murine phenotypic data was available, 37.3% result in pre-weaning lethality (< 3 weeks of age)
length(which(universe_df$mouse_ko=="Y"))

#assessing 15903/19196 (83%) of protein coding genes among eleven different cell lines 
length(which(universe_df$cell_ko=="Y"))/length(universe_df$gene)

#Only 21% of cell essential genes (478/2233) are linked currently to human disorders. 
length(which(universe_df$cell_essential=="Y"&universe_df$omim=="Y"))/length(which(universe_df$cell_essential=="Y"))


#52% of known clinically relevant OMIM genes are classified as tolerant to genetic variation 
length(which(universe_df$omim=="Y"&universe_df$any_constraint=="N"))/length(which(universe_df$omim=="Y"))


#Autosomal recessive (AR) genes rarely demonstrate missense constraint (4.3%) or intolerance to LoF variation (8.4%). 
length(which(universe_df$Inheritance_pattern=="AR"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))
length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))

#genes associated with autosomal dominant (AD) disorders show significantly higher levels of missense constraint (29%), with 46% intolerant to LoF variation

length(which(universe_df$Inheritance_pattern=="AD"&universe_df$mis_z>=3.09))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))


#X-linked (XL) genes show intermediate levels of missense constraint (20.1% XL compared to 4.3% AR and 29% AD); though are the most intolerant of LoF variation (66.9%)
length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&universe_df$mis_z>=3.09))/length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&!is.na(universe_df$exac)))
length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&universe_df$pLI>=0.9))/length(which((grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern))&!is.na(universe_df$exac)))

#59% of OMIM genes linked to early lethality (developmental or neonatal) in knockout mice, compared to 37% non-disease genes
#among our curated list of 624 genes associated with human prenatal, perinatal or infantile lethality – 81% were also associated with a developmental (or neonatal) lethal murine phenotype
length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$mouse_ko=="Y"))
length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$mouse_ko=="Y"))

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

#Unexpectedly, 54 genes classified as cell-essential were not linked to murine lethal phenotypes 
#X/54 genes were associated with sub-viable murine phenotypes and Y/54 with phenotypic abnormality of a certain cell type or organ 
length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"))
xwlength(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"&grepl("cellular phenotype",universe_df$high_MP_phen)))

universe_df$high_MP_phen[which(is.na(universe_df$omim)&universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N")]



