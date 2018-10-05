source("plot_functions.R")

###Figure 1A: Venn diagram showing overlap between cell ‘essentialomes’ ####
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
png('output/Figures/1A.png',width=30,height=30,units="cm",res=1000)
draw.triple.venn(area1 = length(which(common_genes$wang_essentialome=="essential")), area2 = length(which(common_genes$blomen_essentialome=="essential")), area3 = length(which(common_genes$hart_essentialome=="essential")), 
                 n12 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential")),
                 n23 = length(which(common_genes$hart_essentialome=="essential"&common_genes$blomen_essentialome=="essential")), 
                 n13 = length(which(common_genes$wang_essentialome=="essential"&common_genes$hart_essentialome=="essential")), 
                 n123 = length(which(common_genes$wang_essentialome=="essential"&common_genes$blomen_essentialome=="essential"&common_genes$hart_essentialome=="essential")), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("#f1c266", "#f166c5", "#65bbc4"), euler.d = TRUE, scaled = TRUE,fontfamily=rep("Helvetica",7),cat.fontfamily = rep("Helvetica", 3),cex=rep(5,7),cat.cex=rep(3,3))

dev.off()
rm(wang,hart,blomen,common_genes,overrideTriple,a)

###Figure 1B: Number of genes classified as cell-essential among eleven cell lines  ####
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
ggsave("output/Figures/1B.png",height=6.5, width=6.5, units='cm')
dev.off()
rm(cell_line_cutoff_numbers,cell_line_cutoff_levels,cutoffs,cell_cutoff_plot)




###Figure 1C: Venn diagram showing overlap between mouse lethal genes ####
png('output/Figures/1C.png',width=30,height=30,units="cm",res=1000)
draw.pairwise.venn(area1 = length(which(universe_df$lethal_IMPC=="Y")),
                   area2 = length(which(universe_df$lethal_MGI=="Y")), 
                   cross.area=length(which(universe_df$lethal_IMPC=="Y"&universe_df$lethal_MGI=="Y")), 
                   category = c("IMPC", "MGI"), lty = "blank",fill = c("#ef7b0b", "#082e66"),alpha=0.6, euler.d = TRUE, 
                   scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica",cat.cex=(c(5,5)),
                   cat.fontfamily="Helvetica",cat.pos=c(220,140),cat.dist=c(.03,.03))

dev.off()



###Figure 1D: Histograms comparing proportion of protein-coding genes found to be essential in Yeast, Human cells, Mice and humans####
yeast_viable <-read.table("Gene_lists/EG/yeast/viable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable <-read.table("Gene_lists/EG/yeast/inviable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable<-yeast_inviable[which(yeast_inviable$V6=="null"),]
yeast_viable<-yeast_viable[which(yeast_viable$V6=="null"),]

lethal_prop <- data.frame(Yeast=c(length(unique(yeast_inviable$V1))/(length(unique(yeast_inviable$V1))+length(unique(yeast_viable$V1))),
                                  length(unique(yeast_viable$V1))/(length(unique(yeast_inviable$V1))+length(unique(yeast_viable$V1)))),
                          Cells=c(length(which(universe_df$cell_essential_hits>=2))/length(which(!is.na(universe_df$cell_ko))),
                                  length(which(universe_df$cell_essential_hits<2))/length(which(!is.na(universe_df$cell_ko)))),
                          Mouse=c(length(which(universe_df$lethal_mouse=="Y"))/length(which(universe_df$mouse_ko=="Y")),
                                  length(which(universe_df$lethal_mouse=="N"))/length(which(universe_df$mouse_ko=="Y"))),
                          Human=c(length(which(universe_df$human_lethal_B=="Y"))/length(universe_df$gene),
                                  length(which(universe_df$human_lethal_B=="N"))/length(universe_df$gene)))
lethal_propm <- melt(lethal_prop)
lethal_propm$mis_constraint <- rep(c("Lethal     ","Non-lethal     "), 4)
levels <- c(paste0("Yeast \n n = ",(length(unique(yeast_inviable$V1))+length(unique(yeast_viable$V1)))),
            paste0("Human Cells \n n = ",length(which(!is.na(universe_df$cell_ko)))),
            paste0("Mouse \n n = ",length(which(universe_df$mouse_ko=="Y"))),
            paste0("Human \n n = ",length(universe_df$gene)))
lethal_propm$variable <- rep(levels,each=2)
lethal_propm$variable <- factor(lethal_propm$variable, levels = levels)


labels <- rep("",8)
labels[c(1,3,5,7)]<-c(paste0(round(lethal_propm$value[1]*100,0),"%"),
                      paste0(round(lethal_propm$value[3]*100,0),"%"),
                      paste0(round(lethal_propm$value[5]*100,0),"%"),
                      paste0(round(lethal_propm$value[7]*100,0),"%"))
f <- ggplot(dat=lethal_propm, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="black",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("black","steelblue3"))+theme(legend.position="right")+geom_text(aes(label = labels),position = "identity",vjust=-1)
ggsave("output/Figures/1D.pdf",height=9, width=20.5, units='cm')
rm(yeast_inviable,yeast_viable,lethal_prop,lethal_propm,f,labels)






###Figure 1E: Venn diagrams showing overlap between OMIM genes and mouse lethal and cell essential genes####
overrideTriple<- "ye"
#overlap between OMIM, cell + mouse lethality when using cutoff of 3 cell lines- only genes with cell or mouse data
png('output/Figures/1Ei.png',width=30,height=30,units="cm",res=1000)
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
png('output/Figures/1Eii.png',width=10,height=10,units="cm",res=1000)
plot(a,
     quantities = TRUE,
     fill = c("slategray","firebrick4", "sandybrown"),alpha=0.8,
     lty = 0,
     labels=c("OMIM", "Mouse Lethal", "Cell Essential"),fontfamily = "Helvetica")
dev.off()
rm(a,overrideTriple)
