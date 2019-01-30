###Figure 1F: Venn diagrams showing overlap between OMIM genes and mouse lethal and cell essential genes####
overrideTriple<- "ye"
#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- only genes with cell or mouse data
png('Analysis/reviewer_comments/1Fi_1cellline.png',width=,8,height=8,units="cm",res=2000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y"&(universe_df$cell_ko=="Y"|universe_df$mouse_ko=="Y"))),
                 area2 = length(which(universe_df$lethal_mouse=="Y")),
                 area3 = length(which(universe_df$cell_essential_hits>=1)), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=1)), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential_hits>=1)), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=1)), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(2,7),cat.cex=rep(1,3))
dev.off()

png('Analysis/reviewer_comments/1Fi_2celllines.png',width=8,height=8,units="cm",res=2000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y"&(universe_df$cell_ko=="Y"|universe_df$mouse_ko=="Y"))),
                 area2 = length(which(universe_df$lethal_mouse=="Y")),
                 area3 = length(which(universe_df$cell_essential_hits>=2)), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=2)), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential_hits>=2)), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=2)), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(2,7),cat.cex=rep(1,3))
dev.off()
