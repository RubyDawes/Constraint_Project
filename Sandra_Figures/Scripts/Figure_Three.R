
###3a
#get wang, blomen, hart genes and essentialomes
#read wang genes and which are essential in different cell lines
wang <- read.xlsx("Gene_lists/Cell_KOs/wang_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
wang$Gene<- checkGeneSymbols(wang$Gene, unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
wang<- wang[-which(duplicated(wang$Gene)=="TRUE"),] 

#read blomen genes and which are essential in different cell lines
blomen <- read.xlsx("Gene_lists/Cell_KOs/blomen_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
blomen$Gene<- checkGeneSymbols(blomen$Gene, unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
blomen<- blomen[-which(duplicated(blomen$Gene)=="TRUE"),] 

#read hart genes and which are essential in different cell lines
hart <- read.xlsx("Gene_lists/Cell_KOs/hart_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
hart$Gene<- checkGeneSymbols(hart$Gene, unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]
hart<- hart[-which(duplicated(hart$Gene)=="TRUE"),] 

common_genes <-data.frame(intersect(wang$Gene,intersect(hart$Gene,blomen$Gene)))
names(common_genes) <- "genes"
common_genes$wang_essentialome <- vlookup(common_genes$genes,wang,result_column = "essentialome",lookup_column = "Gene")
common_genes$blomen_essentialome <- vlookup(common_genes$genes,blomen,result_column = "essentialome",lookup_column = "Gene")
common_genes$hart_essentialome <- vlookup(common_genes$genes,hart,result_column = "essentialome",lookup_column = "Gene")

overrideTriple<- "ye"
png('Sandra_Figures/Figs/fig3a-BWH_venn_nolabels.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(common_genes$wang_essentialome=="essential")), area2 = length(which(common_genes$blomen_essentialome=="essential")), area3 = length(which(common_genes$hart_essentialome=="essential")), 
                 n12 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential")),
                 n23 = length(which(common_genes$hart_essentialome=="essential"&common_genes$blomen_essentialome=="essential")), 
                 n13 = length(which(common_genes$wang_essentialome=="essential"&common_genes$hart_essentialome=="essential")), 
                 n123 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential"&common_genes$hart_essentialome=="essential")), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("#f1c266", "#f3a292", "#65bbc4"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(0,3))
dev.off()
png('Sandra_Figures/Figs/fig3a-BWH_venn_withlabels.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(common_genes$wang_essentialome=="essential")), area2 = length(which(common_genes$blomen_essentialome=="essential")), area3 = length(which(common_genes$hart_essentialome=="essential")), 
                 n12 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential")),
                 n23 = length(which(common_genes$hart_essentialome=="essential"&common_genes$blomen_essentialome=="essential")), 
                 n13 = length(which(common_genes$wang_essentialome=="essential"&common_genes$hart_essentialome=="essential")), 
                 n123 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential"&common_genes$hart_essentialome=="essential")), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("#f1c266", "#f3a292", "#65bbc4"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))

dev.off()
rm(wang,hart,blomen,common_genes,overrideTriple)


###3b
cell_line_cutoff_levels <- seq(1, 11, by=1)
cell_line_cutoff_numbers <- lapply(cell_line_cutoff_levels,function(x) length(which(universe_df$cell_essential_hits >= x)))
cutoffs <- data.frame(cell_line_cutoff_levels,unlist(cell_line_cutoff_numbers))
names(cutoffs)<-c("levels","numbers")
chosen <- subset(cutoffs, levels== "3")
png('Sandra_Figures/Figs/fig3b_cell_cutoffs.png',width=30,height=30,units="cm",res=1000)
cell_cutoff_plot <- ggplot(cutoffs,aes(x=levels,y=numbers))+
  geom_point(size=6,color="sandybrown")+scatter_theme()+
  labs(x="Number of cell lines",y="Number of essential genes")+
  scale_y_continuous(breaks = pretty(cutoffs$numbers, n = 10))+
  scale_x_continuous(breaks = pretty(cutoffs$levels, n = 10))
ggsave("Sandra_Figures/Figs/fig3b_cell_cutoffs.pdf",height=6.5, width=6.5, units='in')
dev.off()
rm(cell_line_cutoff_numbers,cell_line_cutoff_levels,cutoffs,chosen,cell_cutoff_plot)


###3c
#overlap between cell + mouse lethality when using cutoff of 3 cell lines

mouse=textGrob("Mouse Lethal", gp=gpar(fontsize=40,fontfamily="Helvetica"),rot=90)
cell=textGrob("Cell Essential", gp=gpar(fontsize=40,fontfamily="Helvetica"),rot=270)
a3=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3))+
                        length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)),
                      area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3))+
                        length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), 
                      cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), 
                      category = c("", ""), lty = "blank",fill = c("sandybrown", "firebrick4"), euler.d = TRUE, 
                      scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica")

png('Sandra_Figures/Figs/fig3c_cellmouseoverlap.png',width=30,height=30,units="cm",res=1000)
grid.arrange(gTree(children=a3),left=mouse,right=cell)
dev.off()


