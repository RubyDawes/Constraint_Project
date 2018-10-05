source("plot_functions.R")

#3A: Flow chart linking disease involvement with genetic constraint and murine non-viable phenotypes####
flowchart <- data.frame(OMIM = c("Y","","","","","","","N","","","","","",""),
                        number= c(length(which(universe_df$omim=="Y")),"","","","","","",
                                 length(which(is.na(universe_df$omim))),"","","","","",""),
                        constraint = c("Y","","","N","","","no data","Y","","","N","","","no data"),
                        number = c(length(which(universe_df$omim=="Y"&universe_df$constrained=="Y")),"","",
                                       length(which(universe_df$omim=="Y"&universe_df$constrained=="N")),"","",
                                       length(which(universe_df$omim=="Y"&is.na(universe_df$exac))),
                                       length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),"","",
                                       length(which(is.na(universe_df$omim)&universe_df$constrained=="N")),"","",
                                       length(which(is.na(universe_df$omim)&is.na(universe_df$exac)))),
                        "mouse lethal" = c("Y","N","no data","Y","N","no data","","Y","N","no data",
                                "Y","N","no data",""),
                        number = c(length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y")),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="N")),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&is.na(universe_df$lethal_mouse))),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&universe_df$lethal_mouse=="Y")),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&universe_df$lethal_mouse=="N")),
                                           length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&is.na(universe_df$lethal_mouse))),"",
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&universe_df$lethal_mouse=="N")),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&is.na(universe_df$lethal_mouse))),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&universe_df$lethal_mouse=="Y")),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&universe_df$lethal_mouse=="N")),
                                           length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&is.na(universe_df$lethal_mouse))),""))

g <- tableGrob(flowchart)
ggplot2::ggsave('output/Figures/3A.pdf',g,height=12.5, width=14.5, units='cm')

rm(flowchart,g)

#3B: Odds ratio (OR) analysis linking disease involvement to murine and cellular non-viability ####

OR_test <- function(a,b,c,d){
  matrix <- matrix(c(a,b,c,d),nrow=2,dimnames=list(metric1 = c("Y","N"),metric2 = c("Y","N")))
  test <- fisher.test(matrix,alternative = "two.sided")
  return(test)
}


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
scaleFUN <- function(x) sprintf("%.1f", x)
a<-ggplot(omim_metrics_OR) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill=c("indianred3","black","sandybrown")) +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="grey40", alpha=0.9, size=1)+
  bar_theme()+geom_hline(yintercept=1)+xlab("")+
  scale_y_continuous(trans="log10", limits=c(1,2.6),breaks = trans_breaks('log10', function(x) 10^x),labels = scaleFUN)+
  ylab("Odds Ratio \n (OMIM vs non-OMIM)")+theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank())


#comparing OR, grouped into: AR/MT/XLr AND AD/XLd
omim_constraint_ADXLD_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$constrained=="Y")),
                                   length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                                   length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$constrained=="N")),
                                   length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))

omim_mouseleth_ADXLD_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$lethal_mouse=="Y")),
                                  length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                                  length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$lethal_mouse=="N")),
                                  length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N")))
omim_celless_ADXLD_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$cell_essential=="Y")),
                                length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y")),
                                length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&universe_df$cell_essential=="N")),
                                length(which(is.na(universe_df$omim)&universe_df$cell_essential=="N")))
omim_metrics_ADXLD_OR<-data.frame(categories = c("Constraint","Murine Lethality", "Cell Essentiality"),
                                  OR = c(omim_constraint_ADXLD_OR$estimate,omim_mouseleth_ADXLD_OR$estimate,omim_celless_ADXLD_OR$estimate),
                                  confint_lower = c(omim_constraint_ADXLD_OR$conf.int[1],omim_mouseleth_ADXLD_OR$conf.int[1],omim_celless_ADXLD_OR$conf.int[1]),
                                  confint_upper = c(omim_constraint_ADXLD_OR$conf.int[2],omim_mouseleth_ADXLD_OR$conf.int[2],omim_celless_ADXLD_OR$conf.int[2]))

omim_metrics_ADXLD_OR$categories <- factor(omim_metrics_OR$categories, levels = omim_metrics_OR$categories)
b<-ggplot(omim_metrics_ADXLD_OR) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill=c("indianred3","black","sandybrown")) +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="grey40", alpha=0.9, size=1)+
  bar_theme()+geom_hline(yintercept=1)+xlab("")+
  scale_y_continuous(trans="log10", limits=c(1,5.5),breaks = trans_breaks('log10', function(x) 10^x),labels = scaleFUN)+
  ylab("Odds Ratio \n (AD/XLd OMIM vs non-OMIM)")+theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank())



omim_constraint_ARMTXLR_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&universe_df$constrained=="Y")),
                                     length(which(is.na(universe_df$omim)&universe_df$constrained=="Y")),
                                     length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&universe_df$constrained=="N")),
                                     length(which(is.na(universe_df$omim)&universe_df$constrained=="N")))
omim_mouseleth_ARMTXLR_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&universe_df$lethal_mouse=="Y")),
                                    length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y")),
                                    length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&universe_df$lethal_mouse=="N")),
                                    length(which(is.na(universe_df$omim)&universe_df$lethal_mouse=="N")))
omim_celless_ARMTXLR_OR<- OR_test(length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&universe_df$cell_essential=="Y")),
                                  length(which(is.na(universe_df$omim)&universe_df$cell_essential=="Y")),
                                  length(which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")&universe_df$cell_essential=="N")),
                                  length(which(is.na(universe_df$omim)&universe_df$cell_essential=="N")))
omim_metrics_ARMTXLR_OR<-data.frame(categories = c("Constraint","Murine Lethality", "Cell Essentiality"),
                                    OR = c(omim_constraint_ARMTXLR_OR$estimate,omim_mouseleth_ARMTXLR_OR$estimate,omim_celless_ARMTXLR_OR$estimate),
                                    confint_lower = c(omim_constraint_ARMTXLR_OR$conf.int[1],omim_mouseleth_ARMTXLR_OR$conf.int[1],omim_celless_ARMTXLR_OR$conf.int[1]),
                                    confint_upper = c(omim_constraint_ARMTXLR_OR$conf.int[2],omim_mouseleth_ARMTXLR_OR$conf.int[2],omim_celless_ARMTXLR_OR$conf.int[2]))

omim_metrics_ARMTXLR_OR$categories <- factor(omim_metrics_OR$categories, levels = omim_metrics_OR$categories)

c<-ggplot(omim_metrics_ARMTXLR_OR) +
  geom_bar( aes(x=categories, y=OR), stat="identity", fill=c("indianred3","black","sandybrown")) +
  geom_errorbar( aes(x=categories, ymin=confint_lower, ymax=confint_upper), width=0.2, colour="grey40", alpha=0.9, size=1)+
  bar_theme()+geom_hline(yintercept=1)+xlab("")+
  scale_y_continuous(trans="log10", limits=c(0.55,2.5),breaks = trans_breaks('log10', function(x) 10^x),labels = scaleFUN)+
  ylab("Odds Ratio \n (AR/XLr OMIM vs non-OMIM)")+theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank())


ggarrange(a,b,c,ncol=3)
ggsave("output/Figures/3B.pdf",height=9, width=17, units='cm')

