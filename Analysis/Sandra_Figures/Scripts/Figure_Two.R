###2a
omim_ko <- data.frame(group=c("Available","Unavailable"),
                      value=c(length(which(universe_df$omim=="Y"&!is.na(universe_df$mouse_ko))),
                              length(which(universe_df$omim=="Y"&is.na(universe_df$mouse_ko)))))
nonomim_ko <- data.frame(group=c("Available","Unavailable"),
                      value=c(length(which(is.na(universe_df$omim)&!is.na(universe_df$mouse_ko))),
                              length(which(is.na(universe_df$omim)&is.na(universe_df$mouse_ko)))))
omim_ko$group <- factor(omim_ko$group, levels = omim_ko$group)
nonomim_ko$group <- factor(nonomim_ko$group, levels = nonomim_ko$group)
omim_ko_pie <- ggplot(omim_ko, aes(x="", y=value, fill=group))
omim_ko_pie<- omim_ko_pie+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
omim_ko_pie<- omim_ko_pie+scale_fill_manual(values=c("slategray1","slategray"))+theme(legend.position="none")
omim_ko_pie<- omim_ko_pie+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))
omim_ko_pie<- omim_ko_pie+ggtitle(paste("OMIM genes \n n=",length(which(universe_df$omim=="Y"))))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))

nonomim_ko_pie <- ggplot(nonomim_ko, aes(x="", y=value, fill=group))
nonomim_ko_pie<- nonomim_ko_pie+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
nonomim_ko_pie<- nonomim_ko_pie+scale_fill_manual(values=c("slategray1","slategray"))+theme(legend.position="none")
nonomim_ko_pie<- nonomim_ko_pie+geom_text(aes(y = value,label = percent(value/sum(value))),size=4,position = position_stack(vjust = 0.5))
nonomim_ko_pie<- nonomim_ko_pie+ggtitle(paste("Non-OMIM genes \n n=",length(which(is.na(universe_df$omim)))))+theme(plot.title = element_text(size = 12, face = "bold"),legend.text=element_text(size=12))


ggarrange(omim_ko_pie,nonomim_ko_pie,common.legend=TRUE,legend="bottom")
  ggsave("Analysis/Sandra_Figures/Figs/fig2a-KOpies.pdf",height=6, width=17.8, units='cm')
  dev.off()
rm(omim_ko,nonomim_ko,omim_ko_pie,nonomim_ko_pie)



###2b
#get wang, blomen, hart genes and essentialomes
#read wang genes and which are essential in different cell lines
wang <- read.xlsx("Gene_lists/Cell_KOs/wang_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
wang$Gene<- checkGeneSymbols(wang$Gene, unmapped.as.na=FALSE)[[3]]
wang<- wang[-which(duplicated(wang$Gene)=="TRUE"),] 

#read blomen genes and which are essential in different cell lines
blomen <- read.xlsx("Gene_lists/Cell_KOs/blomen_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
blomen$Gene<- checkGeneSymbols(blomen$Gene, unmapped.as.na=FALSE)[[3]]
blomen<- blomen[-which(duplicated(blomen$Gene)=="TRUE"),] 

#read hart genes and which are essential in different cell lines
hart <- read.xlsx("Gene_lists/Cell_KOs/hart_essentials.xlsx",startRow = 1, colNames = TRUE, rowNames = FALSE)
hart$Gene<- checkGeneSymbols(hart$Gene, unmapped.as.na=FALSE)[[3]]
hart<- hart[-which(duplicated(hart$Gene)=="TRUE"),] 

a<-intersect(wang$Gene,intersect(hart$Gene,blomen$Gene))
common_genes <-data.frame(genes=a[which(a%in%universe_df$gene)])

common_genes$wang_essentialome <- vlookup(common_genes$genes,wang,result_column = "essentialome",lookup_column = "Gene")
common_genes$blomen_essentialome <- vlookup(common_genes$genes,blomen,result_column = "essentialome",lookup_column = "Gene")
common_genes$hart_essentialome <- vlookup(common_genes$genes,hart,result_column = "essentialome",lookup_column = "Gene")

overrideTriple<- "ye"
png('Analysis/Sandra_Figures/Figs/fig3a-BWH_venn_nolabels.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(common_genes$wang_essentialome=="essential")), area2 = length(which(common_genes$blomen_essentialome=="essential")), area3 = length(which(common_genes$hart_essentialome=="essential")), 
                 n12 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential")),
                 n23 = length(which(common_genes$hart_essentialome=="essential"&common_genes$blomen_essentialome=="essential")), 
                 n13 = length(which(common_genes$wang_essentialome=="essential"&common_genes$hart_essentialome=="essential")), 
                 n123 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential"&common_genes$hart_essentialome=="essential")), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("#f1c266", "#f166c5", "#66f1c0"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(0,3))
dev.off()
png('Analysis/Sandra_Figures/Figs/fig3a-BWH_venn_withlabels.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(common_genes$wang_essentialome=="essential")), area2 = length(which(common_genes$blomen_essentialome=="essential")), area3 = length(which(common_genes$hart_essentialome=="essential")), 
                 n12 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential")),
                 n23 = length(which(common_genes$hart_essentialome=="essential"&common_genes$blomen_essentialome=="essential")), 
                 n13 = length(which(common_genes$wang_essentialome=="essential"&common_genes$hart_essentialome=="essential")), 
                 n123 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential"&common_genes$hart_essentialome=="essential")), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("#f1c266", "#f166c5", "#65bbc4"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))

dev.off()
rm(wang,hart,blomen,common_genes,overrideTriple)


###2c
cell_line_cutoff_levels <- seq(1, 11, by=1)
cell_line_cutoff_numbers <- lapply(cell_line_cutoff_levels,function(x) length(which(universe_df$cell_essential_hits >= x)))
cutoffs <- data.frame(cell_line_cutoff_levels,unlist(cell_line_cutoff_numbers))
names(cutoffs)<-c("levels","numbers")
chosen <- subset(cutoffs, levels== "3")

cell_cutoff_plot <- ggplot(cutoffs,aes(x=levels,y=numbers))+
  geom_point(size=3,color="sandybrown")+scatter_theme()+
  labs(x="Number of cell lines",y="Number of essential genes")+
  scale_y_continuous(breaks = pretty(cutoffs$numbers, n = 12),limits=c(0,5000),expand = c(0,0))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11),limits=c(0,11.5),expand = c(0,0))
cell_cutoff_plot
ggsave("Analysis/Sandra_Figures/Figs/fig3b_cell_cutoffs.pdf",height=6.5, width=6.5, units='cm')
dev.off()
rm(cell_line_cutoff_numbers,cell_line_cutoff_levels,cutoffs,chosen,cell_cutoff_plot)


###2d
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

overrideTriple<- "ye"
#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines
png('Analysis/Sandra_Figures/Figs/fig-omim_cell_mouse_venn-LABELS.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y")),
                 area2 = length(which(universe_df$lethal_mouse=="Y")),
                 area3 = length(which(universe_df$cell_essential=="Y")), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))
dev.off()
png('Analysis/Sandra_Figures/Figs/fig-omim_cell_mouse_venn.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y")), area2 = length(which(universe_df$lethal_mouse=="Y")), area3 = length(which(universe_df$cell_essential=="Y")), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(0,3))
dev.off()

#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- only genes with cell or mouse data
png('Analysis/Sandra_Figures/Figs/fig-omim_cell_mouse_venn.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y"&(universe_df$cell_ko=="Y"|universe_df$mouse_ko=="Y"))),
                 area2 = length(which(universe_df$lethal_mouse=="Y")),
                 area3 = length(which(universe_df$cell_essential=="Y")), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))
dev.off()
png('Analysis/Sandra_Figures/Figs/fig-omim_cell_mouse_venn.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y")), area2 = length(which(universe_df$lethal_mouse=="Y")), area3 = length(which(universe_df$cell_essential=="Y")), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(0,3))
dev.off()


#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- but only genes with mouse data
png('Analysis/Sandra_Figures/Figs/fig-omim_cell_mouse_venn_withmousedata-LABELS.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$cell_essential))), 
                 area2 = length(which(universe_df$lethal_mouse=="Y"&!is.na(universe_df$cell_essential))), 
                 area3 = length(which(universe_df$cell_essential=="Y"&!is.na(universe_df$lethal_mouse))), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&!is.na(universe_df$cell_essential))),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y"&!is.na(universe_df$lethal_mouse))), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(1,7),cat.cex=rep(1,3))
dev.off()


#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- but only genes with mouse data
png('Analysis/Sandra_Figures/Figs/fig-omim_cell_mouse_venn_withmousedata.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y"&!is.na(universe_df$lethal_mouse)&!is.na(universe_df$cell_essential))), area2 = length(which(universe_df$lethal_mouse=="Y"&!is.na(universe_df$cell_essential))), area3 = length(which(universe_df$cell_essential=="Y"&!is.na(universe_df$lethal_mouse))), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&!is.na(universe_df$cell_essential))),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y"&!is.na(universe_df$lethal_mouse))), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(0,7),cat.cex=rep(0,3))
dev.off()


