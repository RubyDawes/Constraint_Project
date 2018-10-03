png('Analysis/Sandra_Figures/Figs/venn_1c.png',width=30,height=30,units="cm",res=1000)
draw.pairwise.venn(area1 = length(which(universe_df$lethal_IMPC=="Y")),
                      area2 = length(which(universe_df$lethal_MGI=="Y")), 
                      cross.area=length(which(universe_df$lethal_IMPC=="Y"&universe_df$lethal_MGI=="Y")), 
                      category = c("", ""), lty = "blank",fill = c("#ef7b0b", "#001d49"), euler.d = TRUE, 
                      scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica")
dev.off()


draw.pairwise.venn(area1 = length(which(!is.na(universe_df$lethal_IMPC))),
                      area2 = length(which(!is.na(universe_df$lethal_MGI))), 
                      cross.area=length(which(!is.na(universe_df$lethal_IMPC)&!is.na(universe_df$lethal_MGI))), 
                      category = c("", ""), lty = "blank",fill = c("midnightblue", "orangered4"), euler.d = TRUE, 
                      scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica")

