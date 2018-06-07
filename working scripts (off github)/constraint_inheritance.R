
source("2.ExAC_constraint/exac_constraint.R")
source("plot_functions.R")

omim$mis_z <- vlookup(omim$gene, exac,result_column="mis_z",lookup_column="gene")
omim$pLI <- vlookup(omim$gene, exac,result_column="pLI",lookup_column="gene")

omim_AD_mis <- omim$mis_z[which(grepl(pattern="AD",x=omim$Inheritance_pattern,ignore.case = FALSE))]
omim_AD_mis <- omim_AD_mis[!is.na(omim_AD_mis)]
omim_AR_mis <- omim$mis_z[which(grepl(pattern="AR",x=omim$Inheritance_pattern,ignore.case = FALSE))]
omim_AR_mis <- omim_AR_mis[!is.na(omim_AR_mis)]
omim_XL_mis <- omim$mis_z[which(grepl(pattern="XL",x=omim$Inheritance_pattern,ignore.case = FALSE))]
omim_XL_mis <- omim_XL_mis[!is.na(omim_XL_mis)]

AD_constraint <- data.frame(category = c("constrained","not constrained"),values = c(length(which(omim_AD_mis>=3.09|omim_AD_lof>=0.9)),length(which(omim_AD_mis<3.09&omim_AD_lof<0.9))))
AD_constraint$percents <- percent(AD_constraint$values/sum(AD_constraint$values))
AR_constraint <- data.frame(category = c("constrained","not constrained"),values = c(length(which(omim_AR_mis>=3.09|omim_AR_lof>=0.9)),length(which(omim_AR_mis<3.09&omim_AR_lof<0.9))))
AR_constraint$percents <- percent(AR_constraint$values/sum(AR_constraint$values))
XL_constraint <- data.frame(category = c("constrained","not constrained"),values = c(length(which(omim_XL_mis>=3.09|omim_XL_lof>=0.9)),length(which(omim_XL_mis<3.09&omim_XL_lof<0.9))))
XL_constraint$percents <- percent(XL_constraint$values/sum(XL_constraint$values))

misdata <- data.frame(category = rep(c("AD","AR","XL"), c(length(omim_AD_mis),length(omim_AR_mis),length(omim_XL_mis))), missense_Z_score = c(omim_AD_mis,omim_AR_mis,omim_XL_mis))
missense<- ggplot(misdata, aes(x=missense_Z_score,color=category,fill=category))
missense <- missense+density_theme()+geom_density(alpha=0.5)
missense <- missense+geom_segment(aes(x=3.09,y=0,xend=3.09,yend=0.3),color="slategray",linetype="dashed",size=0.5)
#missense <- missense +geom_segment(aes(x = mean(omim_AD_mis), y = 0, xend = mean(omim_AD_mis), yend = 0.2),color="dodgerblue3",linetype="dashed",size=0.5,alpha=0.5)
#missense <- missense +geom_segment(aes(x = mean(omim_AR_mis), y = 0, xend = mean(omim_AR_mis), yend = 0.28),color="brown3",linetype="dashed",size=0.5,alpha=0.5)
#missense <- missense +geom_segment(aes(x = mean(omim_XL_mis), y = 0, xend = mean(omim_XL_mis), yend = 0.25),color="snow2",linetype="dashed",size=0.5,alpha=0.5)
missense <- missense+scale_fill_manual(values=c("dodgerblue3","brown3","snow"))+scale_color_manual(values=c("slategray","slategray","slategray"))+coord_cartesian(xlim = c(-5, 10)) 

omim_AD_lof <- omim$pLI[which(grepl(pattern="AD",x=omim$Inheritance_pattern,ignore.case = FALSE))]
omim_AD_lof <- omim_AD_lof[!is.na(omim_AD_lof)]
omim_AR_lof <- omim$pLI[which(grepl(pattern="AR",x=omim$Inheritance_pattern,ignore.case = FALSE))]
omim_AR_lof <- omim_AR_lof[!is.na(omim_AR_lof)]
omim_XL_lof <- omim$pLI[which(grepl(pattern="XL",x=omim$Inheritance_pattern,ignore.case = FALSE))]
omim_XL_lof <- omim_XL_lof[!is.na(omim_XL_lof)]

lofdata <- data.frame(category = rep(c("AD","AR","XL"), c(length(omim_AD_lof),length(omim_AR_lof),length(omim_XL_lof))), LoF_pLI_score = c(omim_AD_lof,omim_AR_lof,omim_XL_lof))
lof<- ggplot(lofdata, aes(x=LoF_pLI_score,color=category,fill=category))+geom_density(alpha=0.2)
lof <- lof+density_theme()+geom_density(alpha=0.5)
lof <- lof+geom_segment(aes(x=0.9,y=0,xend=0.9,yend=11),color="slategray",linetype="dashed",size=0.5)
lof <- lof+scale_fill_manual(values=c("dodgerblue3","brown3","snow"))+scale_color_manual(values=c("slategray","slategray","slategray"))+coord_cartesian(xlim = c(0, 1)) 

png("output/figures/fig4.png",width=1000,height=1000,type="quartz",res=150)
grid_arrange_shared_legend(missense,lof)
dev.off()

