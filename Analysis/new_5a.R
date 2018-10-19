#getting annos
goslim_anno <- read.table("Gene_lists/GOSLIM/annotations.mapped.gaf",header=FALSE,sep="\t",comment.char="!",fill=TRUE)
goslim_anno <- goslim_anno[-which(grepl("NOT",goslim_anno$V4)),]

goslim <- get_ontology("Gene_lists/GOSLIM/goslim_generic.obo", propagate_relationships = "is_a",extract_tags = "minimal")

category <- unlist(lapply(goslim$ancestors[1:149],function(x) ifelse("GO:0003674"%in%x,"MF",
                                                                     ifelse("GO:0008150"%in%x,"BP",
                                                                            ifelse("GO:0005575"%in%x,"CC",NA)))))

go_fish <- data.frame(go_id = goslim$id[1:149],go_name = goslim$name[1:149],category = category)
rownames(go_fish) <- c()

go_fish$annotated_genes <- lapply(go_fish$go_id,function(x) 
  unique(goslim_anno$V3[which(as.character(goslim_anno$V3)%in%universe_df$gene&as.character(goslim_anno$V5)==x)]))
go_fish$gene_no <- lengths(go_fish$annotated_genes)


go_fish$lethal_annotated_genes <-lapply(go_fish$go_id,function(x) 
  unique(goslim_anno$V3[which(as.character(goslim_anno$V3)%in%universe_df$gene[which(universe_df$human_lethal_B=="Y")]&as.character(goslim_anno$V5)==x)]))
go_fish$lethal_gene_no<-lengths(go_fish$lethal_annotated_genes)


#go_fish putative lethal dom vs putative lethal rec####
go_fish$putlethaldom_annotated_genes <- lapply(go_fish$go_id,function(x) 
  unique(goslim_anno$V3[which(as.character(goslim_anno$V3)%in%universe_df$gene[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y")]&as.character(goslim_anno$V5)==x)]))
go_fish$putlethaldom_gene_no <- lengths(go_fish$putlethaldom_annotated_genes)
go_fish$putlethaldom_gene_perc <- go_fish$putlethaldom_gene_no/length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno$V3))

go_fish$putlethalrec_annotated_genes <- lapply(go_fish$go_id,function(x) 
  unique(goslim_anno$V3[which(as.character(goslim_anno$V3)%in%universe_df$gene[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N")]&as.character(goslim_anno$V5)==x)]))
go_fish$putlethalrec_gene_no <- lengths(go_fish$putlethalrec_annotated_genes)
go_fish$putlethalrec_gene_perc <- go_fish$putlethalrec_gene_no/length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno$V3))

#put lethal dom vs rec####
annotated_putdomlethal_no <-length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno$V3))
annotated_putreclethal_no <-length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno$V3))

go_fish$putdomvsreclethal_odds_ratio <- unlist(lapply(1:length(go_fish$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish$putlethaldom_gene_no[x],go_fish$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$estimate[[1]])
}))

go_fish$putdomvsreclethal_confint_lower <- unlist(lapply(1:length(go_fish$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish$putlethaldom_gene_no[x],go_fish$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[1])
}))

go_fish$putdomvsreclethal_confint_higher <- unlist(lapply(1:length(go_fish$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish$putlethaldom_gene_no[x],go_fish$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$conf.int[2])
}))

go_fish$putdomvsreclethal_pval <- unlist(lapply(1:length(go_fish$go_name),function(x) {
  fish<-fisher.test(matrix(c(go_fish$putlethaldom_gene_no[x],go_fish$putlethalrec_gene_no[x],
                             annotated_putdomlethal_no-go_fish$putlethaldom_gene_no[x],annotated_putreclethal_no-go_fish$putlethalrec_gene_no[x]
  ),ncol=2),alternative="two.sided")
  return(fish$p.value)
}))


go_fish <- go_fish[which(go_fish$putdomvsreclethal_pval<0.05),]
go_fish <- go_fish[-which(go_fish$putlethaldom_gene_no<20&go_fish$putlethalrec_gene_no<20),]

#plotting heatmaps####
go_fish$gene <- rep("gene",length(go_fish$go_id))

#putative dom vs rec
go_fish_putlethal_CC <- go_fish[which(go_fish$category=="CC"),]
go_fish_putlethal_CC <- go_fish_putlethal_CC[order(go_fish_putlethal_CC$putdomvsreclethal_odds_ratio),]
go_fish_putlethal_CC$go_name <- factor(go_fish_putlethal_CC$go_name,levels=go_fish_putlethal_CC$go_name)
g<-ggplot(data = go_fish_putlethal_CC, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = putdomvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "colourbar", limits=c(0.1,3))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position = 'none',axis.title.y=element_blank(),axis.text.y=element_text(size=9))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_putlethal_MF <- go_fish[which(go_fish$category=="MF"),]
go_fish_putlethal_MF <- go_fish_putlethal_MF[order(go_fish_putlethal_MF$putdomvsreclethal_odds_ratio),]
go_fish_putlethal_MF$go_name <- factor(go_fish_putlethal_MF$go_name,levels=go_fish_putlethal_MF$go_name)
h<-ggplot(data = go_fish_putlethal_MF, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = putdomvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "colourbar", limits=c(0.1,4))+
  bar_theme()+theme(axis.text.x=element_blank(),legend.position = 'none',axis.title.y=element_blank(),axis.text.y=element_text(size=9))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))

go_fish_putlethal_BP <- go_fish[which(go_fish$category=="BP"),]
go_fish_putlethal_BP <- go_fish_putlethal_BP[order(go_fish_putlethal_BP$putdomvsreclethal_odds_ratio),]
go_fish_putlethal_BP$go_name <- factor(go_fish_putlethal_BP$go_name,levels=go_fish_putlethal_BP$go_name)
i<-ggplot(data = go_fish_putlethal_BP, aes(x= gene,y = go_name)) +
  geom_tile(aes(fill = putdomvsreclethal_odds_ratio,width=1)) +
  scale_fill_gradientn(colours = c("#a020f0","white","#008b45"), 
                       values = rescale(c(0.1,1,3)),breaks=c(0.1,0.5,1.5,2,3),labels=c(0.1,0.5,1.5,2,3),
                       guide = "colourbar", limits=c(0.1,4))+
  bar_theme()+theme(axis.text.x=element_blank(),axis.title.y=element_blank(),axis.text.y=element_text(size=9))+scale_y_discrete(expand = c(0, 0),labels = function(x) lapply(strwrap(x, width = 1, simplify = FALSE), paste, collapse="\n"))+scale_x_discrete(expand = c(0, 0))



