OR_test <- function(a,b,c,d){
  matrix <- matrix(c(a,b,c,d),nrow=2,dimnames=list(metric1 = c("Y","N"),metric2 = c("Y","N")))
  test <- fisher.test(matrix,alternative = "two.sided")
  return(test)
}

#comparing ORs for constraint different inheritances

omim_constraint_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$constrained=="Y")),
                             length(which(universe_df$omim=="Y"&universe_df$constrained=="N")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_MT_constraint_OR<-OR_test(length(which(universe_df$omim=="Y"&grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$constrained=="Y")),
                               length(which(universe_df$omim=="Y"&grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$constrained=="N")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_AR_constraint_OR<-OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="Y")),
                               length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="N")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_ARAD_constraint_OR<-OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR,AD"&universe_df$constrained=="Y")),
                                 length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR,AD"&universe_df$constrained=="N")),
                                 length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                                 length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_AD_constraint_OR<-OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AD"&universe_df$constrained=="Y")),
                               length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AD"&universe_df$constrained=="N")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_XL_constraint_OR<-OR_test(length(which(universe_df$omim=="Y"&grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$constrained=="Y")),
                               length(which(universe_df$omim=="Y"&grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$constrained=="N")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                               length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))

omim_constraint_inheritance_OR<-data.frame(categories = c("OMIM","MT", "AR", "AR,AD","AD","XL"),
                                           OR = c(omim_constraint_OR$estimate,omim_MT_constraint_OR$estimate,omim_AR_constraint_OR$estimate,omim_ARAD_constraint_OR$estimate,omim_AD_constraint_OR$estimate,omim_XL_constraint_OR$estimate),
                                           confint_lower = c(omim_constraint_OR$conf.int[1],omim_MT_constraint_OR$conf.int[1],omim_AR_constraint_OR$conf.int[1],omim_ARAD_constraint_OR$conf.int[1],omim_AD_constraint_OR$conf.int[1],omim_XL_constraint_OR$conf.int[1]),
                                           confint_upper = c(omim_constraint_OR$conf.int[2],omim_MT_constraint_OR$conf.int[2],omim_AR_constraint_OR$conf.int[2],omim_ARAD_constraint_OR$conf.int[2],omim_AD_constraint_OR$conf.int[2],omim_XL_constraint_OR$conf.int[2]))
rm(omim_constraint_OR,omim_MT_constraint_OR,omim_AD_constraint_OR,omim_ARAD_constraint_OR,omim_AR_constraint_OR,omim_XL_constraint_OR)
omim_constraint_inheritance_OR$categories <- factor(omim_constraint_inheritance_OR$categories, levels = omim_constraint_inheritance_OR$categories)

ggplot(omim_constraint_inheritance_OR[-2,]) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill="grey") +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="orange", alpha=0.9, size=1)+
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10", limits=c(NA,13),breaks = c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10))+ylab("Odds Ratio \n (Constrained vs Not constrained)")
ggsave("Analysis/Sandra_Figures/Figs/OR_inheritance.pdf",height=7, width=7, units='cm')

#comparing ORs for mouse lethal, constraint, cell essential
omim_constraint_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$constrained=="Y")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                             length(which(universe_df$omim=="Y"&universe_df$constrained=="N")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_mouseleth_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                            length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                             length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="N")),
                             length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N")))
omim_celless_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")),
                            length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y")),
                             length(which(universe_df$omim=="Y"&universe_df$cell_essential=="N")),
                             length(which(is.na(universe_df$omim)&universe_df$cell_essential=="N")))
omim_metrics_OR<-data.frame(categories = c("Constraint","Murine Lethality", "Cell Essentiality"),
                            OR = c(omim_constraint_OR$estimate,omim_mouseleth_OR$estimate,omim_celless_OR$estimate),
                            confint_lower = c(omim_constraint_OR$conf.int[1],omim_mouseleth_OR$conf.int[1],omim_celless_OR$conf.int[1]),
                            confint_upper = c(omim_constraint_OR$conf.int[2],omim_mouseleth_OR$conf.int[2],omim_celless_OR$conf.int[2]))

omim_metrics_OR$categories <- factor(omim_metrics_OR$categories, levels = omim_metrics_OR$categories)

ggplot(omim_metrics_OR) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill="grey") +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="orange", alpha=0.9, size=1)+
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10", limits=c(NA,3),breaks = c(1.0,1.5, 2.0,2.5,3.0))+
  ylab("Odds Ratio \n (OMIM vs non-OMIM)")+theme(axis.text.x = element_text(angle=60, hjust=1))
ggsave("Analysis/Sandra_Figures/Figs/OR_metrics.pdf",height=7, width=7, units='cm')

#comparing OR for just AD genes
omim_constraint_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AD"&universe_df$constrained=="Y")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                             length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AD"&universe_df$constrained=="N")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_mouseleth_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AD"&universe_df$lethal_mouse=="Y")),
                            length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                            length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AD"&universe_df$lethal_mouse=="N")),
                            length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N")))
omim_celless_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AD"&universe_df$cell_essential=="Y")),
                          length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y")),
                          length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AD"&universe_df$cell_essential=="N")),
                          length(which(is.na(universe_df$omim)&universe_df$cell_essential=="N")))
omim_metrics_AD_OR<-data.frame(categories = c("Constraint","Murine Lethality", "Cell Essentiality"),
                            OR = c(omim_constraint_OR$estimate,omim_mouseleth_OR$estimate,omim_celless_OR$estimate),
                            confint_lower = c(omim_constraint_OR$conf.int[1],omim_mouseleth_OR$conf.int[1],omim_celless_OR$conf.int[1]),
                            confint_upper = c(omim_constraint_OR$conf.int[2],omim_mouseleth_OR$conf.int[2],omim_celless_OR$conf.int[2]))

omim_metrics_AD_OR$categories <- factor(omim_metrics_OR$categories, levels = omim_metrics_OR$categories)

ggplot(omim_metrics_AD_OR) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill="grey") +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="orange", alpha=0.9, size=1)+
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10", limits=c(NA,5.5),breaks = c(1.0,1.5, 2.0,2.5,3.0,4,5))+
  ylab("Odds Ratio \n (AD OMIM vs non-OMIM)")+theme(axis.text.x = element_text(angle=60, hjust=1))
ggsave("Analysis/Sandra_Figures/Figs/OR_AD_metrics.pdf",height=7, width=7, units='cm')


#comparing OR for just AR genes
omim_constraint_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AR"&universe_df$constrained=="Y")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                             length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AR"&universe_df$constrained=="N")),
                             length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_mouseleth_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AR"&universe_df$lethal_mouse=="Y")),
                            length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                            length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AR"&universe_df$lethal_mouse=="N")),
                            length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N")))
omim_celless_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AR"&universe_df$cell_essential=="Y")),
                          length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y")),
                          length(which(universe_df$omim=="Y"&universe_df$Inheritance=="AR"&universe_df$cell_essential=="N")),
                          length(which(is.na(universe_df$omim)&universe_df$cell_essential=="N")))
omim_metrics_AR_OR<-data.frame(categories = c("Constraint","Murine Lethality", "Cell Essentiality"),
                               OR = c(omim_constraint_OR$estimate,omim_mouseleth_OR$estimate,omim_celless_OR$estimate),
                               confint_lower = c(omim_constraint_OR$conf.int[1],omim_mouseleth_OR$conf.int[1],omim_celless_OR$conf.int[1]),
                               confint_upper = c(omim_constraint_OR$conf.int[2],omim_mouseleth_OR$conf.int[2],omim_celless_OR$conf.int[2]))

omim_metrics_AR_OR$categories <- factor(omim_metrics_OR$categories, levels = omim_metrics_OR$categories)

ggplot(omim_metrics_AR_OR) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill="grey") +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="orange", alpha=0.9, size=1)+
  bar_theme_or()+geom_hline(yintercept=1)+
  scale_y_continuous(trans="log10", limits=c(NA,2.5),breaks = c(0.1,0.5,0.7,1.0,1.5, 2.0,2.5,3.0,4,5))+
  ylab("Odds Ratio \n (AR OMIM vs non-OMIM)")+theme(axis.text.x = element_text(angle=60, hjust=1))
ggsave("Analysis/Sandra_Figures/Figs/OR_AR_metrics.pdf",height=7, width=7, units='cm')
