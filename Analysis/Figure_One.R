source("plot_functions.R")

###Figure 1B: Pies comparing proportion of protein-coding genes found to be essential in Yeast, Human cells, Mice and humans####
yeast_viable <-read.table("Gene_lists/EG/yeast/viable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable <-read.table("Gene_lists/EG/yeast/inviable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable<-yeast_inviable[which(yeast_inviable$V6=="null"),]
yeast_viable<-yeast_viable[which(yeast_viable$V6=="null"),]

yeast_lethal_pie <- data.frame(group=c("lethal", "non-lethal"),
                               value=c(
                                 length(unique(yeast_inviable$V1)),
                                 length(unique(yeast_viable$V1))))
yeast_lethal_pie$group <- factor(yeast_lethal_pie$group, levels = yeast_lethal_pie$group)
a <- ggplot(yeast_lethal_pie, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)+
  scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="none",plot.title=element_text(size=10,face="plain"))+
  geom_text(aes(y = value,label = percent(value/sum(value)),colour=group),size=3,position = position_stack(vjust = 0.5),show.legend=FALSE,fontface="bold")+
  scale_colour_manual(values=c("white","black"))+ggtitle(paste0("Yeast \n n = ",sum(yeast_lethal_pie[,2])))
cell_lethal_pie <- data.frame(group=c("lethal", "non-lethal"),
                              value=c(
                                length(which(universe_df$cell_essential_hits>=2)),
                                length(which(universe_df$cell_essential_hits<2))))
cell_lethal_pie$group <- factor(cell_lethal_pie$group, levels = cell_lethal_pie$group)
b <- ggplot(cell_lethal_pie, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)+
  scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="none",plot.title=element_text(size=10,face="plain"))+
  geom_text(aes(y = value,label = percent(value/sum(value)),colour=group),size=3,position = position_stack(vjust = 0.5),show.legend=FALSE,fontface="bold")+
  scale_colour_manual(values=c("white","black"))+ggtitle(paste0("Human Cells \n n = ",sum(cell_lethal_pie[,2])))
mouse_lethal_pie <- data.frame(group=c("lethal", "non-lethal"),
                              value=c(
                                length(which(universe_df$lethal_mouse=="Y")),
                                length(which(universe_df$lethal_mouse=="N"))))
mouse_lethal_pie$group <- factor(mouse_lethal_pie$group, levels = mouse_lethal_pie$group)
c <- ggplot(mouse_lethal_pie, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)+
  scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="none",plot.title=element_text(size=10,face="plain"))+
  geom_text(aes(y = value,label = percent(value/sum(value)),colour=group),size=3,position = position_stack(vjust = 0.5),show.legend=FALSE,fontface="bold")+
  scale_colour_manual(values=c("white","black"))+ggtitle(paste0("Mouse \n n = ",sum(mouse_lethal_pie[,2])))
human_lethal_pie <- data.frame(group=c("lethal", "non-lethal"),
                               value=c(
                                 length(which(universe_df$human_lethal_B=="Y")),
                                 length(which(universe_df$human_lethal_B=="N"))))
human_lethal_pie$group <- factor(human_lethal_pie$group, levels = human_lethal_pie$group)
d <- ggplot(human_lethal_pie, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)+
  scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="none",plot.title=element_text(size=10,face="plain"))+
  geom_text(aes(y = value,label = percent(value/sum(value)),colour=group),size=3,position = position_stack(vjust = 0.5),show.legend=FALSE,fontface="bold")+
  scale_colour_manual(values=c("white","black"))+ggtitle(paste0("Human \n n = ",sum(human_lethal_pie[,2])))

ggarrange(a,b,c,d,ncol=4,nrow=1)
ggsave("output/Figures/1B.pdf",height=4.5, width=18, units='cm')
rm(a,b,c,d,yeast_inviable,yeast_viable,yeast_lethal_pie,cell_lethal_pie,mouse_lethal_pie,human_lethal_pie)

###Figure 1C: Venn diagram showing overlap between mouse lethal genes ####
png('output/Figures/1C.png',width=30,height=30,units="cm",res=1000)
draw.pairwise.venn(area1 = length(which(universe_df$lethal_IMPC=="Y")),
                   area2 = length(which(universe_df$lethal_MGI=="Y")), 
                   cross.area=length(which(universe_df$lethal_IMPC=="Y"&universe_df$lethal_MGI=="Y")), 
                   category = c("IMPC", "MGI"), lty = "blank",fill = c("#ef7b0b", "#082e66"),alpha=0.6, euler.d = TRUE, 
                   scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica",cat.cex=(c(5,5)),
                   cat.fontfamily="Helvetica",cat.pos=c(220,140),cat.dist=c(.03,.03))

dev.off()




###Figure 1D: Venn diagram showing overlap between cell ‘essentialomes’ ####
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
png('output/Figures/1D.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(common_genes$wang_essentialome=="essential")), area2 = length(which(common_genes$blomen_essentialome=="essential")), area3 = length(which(common_genes$hart_essentialome=="essential")), 
                 n12 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential")),
                 n23 = length(which(common_genes$hart_essentialome=="essential"&common_genes$blomen_essentialome=="essential")), 
                 n13 = length(which(common_genes$wang_essentialome=="essential"&common_genes$hart_essentialome=="essential")), 
                 n123 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential"&common_genes$hart_essentialome=="essential")), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("#f1c266", "#f166c5", "#65bbc4"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))

dev.off()
rm(wang,hart,blomen,common_genes,overrideTriple,a)



###Figure 1E: Number of genes classified as cell-essential among eleven cell lines  ####
cell_line_cutoff_levels <- seq(1, 11, by=1)
cell_line_cutoff_numbers <- lapply(cell_line_cutoff_levels,function(x) length(which(universe_df$cell_essential_hits >= x)))
cutoffs <- data.frame(cell_line_cutoff_levels,unlist(cell_line_cutoff_numbers))
names(cutoffs)<-c("levels","numbers")

cell_cutoff_plot <- ggplot(cutoffs,aes(x=levels,y=numbers))+
  geom_point(size=3,color="sandybrown")+scatter_theme()+
  labs(x="Number of cell lines",y="Number of essential genes")+
  scale_y_continuous(breaks = pretty(cutoffs$numbers, n = 12),limits=c(0,5000),expand = c(0,0))+
  scale_x_continuous(breaks = c(0,1,2,3,4,5,6,7,8,9,10,11),limits=c(0,11.5),expand = c(0,0))
cell_cutoff_plot
ggsave("output/Figures/1E.png",height=6.5, width=6.5, units='cm')
dev.off()
rm(cell_line_cutoff_numbers,cell_line_cutoff_levels,cutoffs,cell_cutoff_plot)




###Figure 1F: Venn diagrams showing overlap between OMIM genes and mouse lethal and cell essential genes####
overrideTriple<- "ye"
#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- only genes with cell or mouse data
png('output/Figures/1Fi.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(universe_df$omim=="Y"&(universe_df$cell_ko=="Y"|universe_df$mouse_ko=="Y"))),
                 area2 = length(which(universe_df$lethal_mouse=="Y")),
                 area3 = length(which(universe_df$cell_essential=="Y")), 
                 n12 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y")),
                 n23 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), 
                 n13 = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y")), 
                 n123 = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")), category = c("OMIM", "Mouse Lethal", "Cell Essential"), lty = "blank", 
                 fill = c("slategray","firebrick4", "sandybrown"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))
dev.off()

#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- but only genes with mouse data- need ellipses to represent areas accurately
library(eulerr)
a<-euler(c("A" = length(which(universe_df$omim=="Y"&universe_df$cell_essential=="N"&universe_df$lethal_mouse=="N")),
        "B" = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N"&is.na(universe_df$omim))), 
        "C" = length(which(universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"&is.na(universe_df$omim))),
        "A&B" = length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")),
        "B&C" = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y"&is.na(universe_df$omim))),
        "A&C" =length(which(universe_df$omim=="Y"&universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N")),
        "A&B&C" =length(which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y"))),shape="ellipse")
png('output/Figures/1Fii.png',width=10,height=10,units="cm",res=1000)
plot(a,
     quantities = TRUE,
     fill = c("slategray","firebrick4", "sandybrown"),alpha=0.8,
     lty = 0,
     labels=c("OMIM", "Mouse Lethal", "Cell Essential"),fontfamily = "Helvetica")
dev.off()
rm(a,overrideTriple)
