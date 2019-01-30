

#are all the genes lethal in heterozygous KO mice, also lethal in homozygous KO mice?
lethal_homhet <- universe_df[which(universe_df$lethal_mouse=="Y"&universe_df$lethal_het_mouse=="Y"),]
lethal_hetonly <- universe_df[which(universe_df$lethal_mouse=="N"&universe_df$lethal_het_mouse=="Y"),]
lethal_hetnohom <- universe_df[which(is.na(universe_df$lethal_mouse)&universe_df$lethal_het_mouse=="Y"),]
lethal_homonly <- universe_df[which(universe_df$lethal_mouse=="Y"&universe_df$lethal_het_mouse=="N"),]
lethal_homnohet <- universe_df[which(universe_df$lethal_mouse=="Y"&is.na(universe_df$lethal_het_mouse)),]

#making venn diagram
###Figure 1C: Venn diagram showing overlap between mouse lethal genes ####
png('Analysis/reviewer comments/mouse_heterozygousKO_lethality_overlap.png',width=30,height=30,units="cm",res=1000)
draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="Y"&!is.na(universe_df$lethal_het_mouse))),
                   area2 = length(which(!is.na(universe_df$lethal_mouse)&universe_df$lethal_het_mouse=="Y")), 
                   cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$lethal_het_mouse=="Y")), 
                   category = c("homozygous KO mice", "heterozygous KO mice"), lty = "blank",fill = c("#ef7b0b", "#082e66"),alpha=0.6, euler.d = TRUE, 
                   scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica",cat.cex=(c(2,2)),
                   cat.fontfamily="Helvetica",cat.pos=c(220,140),cat.dist=c(.03,.03))

dev.off()


length(which(universe_df$mouse_ko=="Y"|universe_df$mouse_het_ko=="Y"))

length(which(is.na(universe_df$mouse_ko)&universe_df$mouse_het_ko=="Y"&universe_df$lethal_het_mouse==))
length(which(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y"))
