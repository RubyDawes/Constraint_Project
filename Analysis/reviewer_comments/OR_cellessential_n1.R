
OR_test <- function(a,b,c,d){
  matrix <- matrix(c(a,b,c,d),nrow=2,dimnames=list(metric1 = c("Y","N"),metric2 = c("Y","N")))
  test <- fisher.test(matrix,alternative = "two.sided")
  return(test)
}

omim_celless_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$cell_essential_hits>=1)),
                          length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits>=1)),
                          length(which(universe_df$omim=="Y"&universe_df$cell_essential_hits==0)),
                          length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits==0)))


omim_celless_ADXLD_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$cell_essential_hits>=1)),
                                length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits>=1)),
                                length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$cell_essential_hits==0)),
                                length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits==0)))


omim_celless_ARMTXLR_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&
                                                 universe_df$cell_essential_hits>=1)),
                                  length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits>=1)),
                                  length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&
                                                 universe_df$cell_essential_hits==0)),
                                  length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits==0)))

lethal_celless_OR<- OR_test(length(which(universe_df$human_lethal_B=="Y"&universe_df$cell_essential_hits>=1)),
                            length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits>=1)),
                            length(which(universe_df$human_lethal_B=="Y"&universe_df$cell_essential_hits==0)),
                            length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits==0)))


lethal_celless_ADXLD_OR<- OR_test(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$cell_essential_hits>=1)),
                                  length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits>=1)),
                                  length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$cell_essential_hits==0)),
                                  length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits==0)))


lethal_celless_ARMTXLR_OR<- OR_test(length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")&universe_df$cell_essential_hits>=1)),
                                    length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits>=1)),
                                    length(which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")&universe_df$cell_essential_hits==0)),
                                    length(which(is.na(universe_df$omim)&universe_df$cell_essential_hits==0)))
