###Figure 1C: Venn diagram showing overlap between mouse lethal genes ####
png('Analysis/reviewer_comments/Figures/1C.png',width=30,height=30,units="cm",res=1000)
draw.pairwise.venn(area1 = length(which(universe_df$lethal_IMPC=="Y"&!is.na(universe_df$lethal_MGI))),
                   area2 = length(which(universe_df$lethal_MGI=="Y"&!is.na(universe_df$lethal_IMPC))), 
                   cross.area=length(which(universe_df$lethal_IMPC=="Y"&universe_df$lethal_MGI=="Y")), 
                   category = c("IMPC", "MGI"), lty = "blank",fill = c("#ef7b0b", "#082e66"),alpha=0.6, euler.d = TRUE, 
                   scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica",cat.cex=(c(5,5)),
                   cat.fontfamily="Helvetica",cat.pos=c(220,140),cat.dist=c(.03,.03))

dev.off()

length(which(universe_df$lethal_IMPC=="N"&!is.na(universe_df$lethal_MGI)))
length(which(universe_df$lethal_IMPC=="N"&universe_df$lethal_MGI=="N"))
length(which(universe_df$lethal_IMPC=="Y"&universe_df$lethal_MGI=="Y"))/
  length(which((universe_df$lethal_IMPC=="Y"|universe_df$lethal_MGI=="Y")&!is.na(universe_df$lethal_IMPC)&!is.na(universe_df$lethal_MGI)))
