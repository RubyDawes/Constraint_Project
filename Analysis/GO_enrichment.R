#getting V3 annos ####
goslim_anno_V3 <- read.table("Gene_lists/GOSLIM/annotations.mapped.V3.gaf",header=FALSE,sep="\t",comment.char="!",fill=TRUE)
goslim_anno_V3$V3 <- checkGeneSymbols(goslim_anno_V3$V3,unmapped.as.na=FALSE)[[3]]
goslim_anno_V3 <- goslim_anno_V3[-which(grepl("NOT",goslim_anno_V3$V4)),]

#go_fish_V3 all vs omim####
go_fish_V3 <- read.table("Gene_lists/GOSLIM/Custom_GOSLIM_V3.txt",header=TRUE,fill=TRUE,sep="\t")
go_fish_V3$go_name<- factor(go_fish_V3$go_name,levels=go_fish_V3$go_name)

go_fish_V3$nonomim_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(is.na(universe_df$omim))]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$gene_no <- lengths(go_fish_V3$nonomim_annotated_genes)
go_fish_V3$gene_perc <- go_fish_V3$gene_no/length(which(universe_df$gene[which(is.na(universe_df$omim))]%in%goslim_anno_V3$V3))


go_fish_V3$omim_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(universe_df$omim=="Y")]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$omim_gene_no <- lengths(go_fish_V3$omim_annotated_genes)
go_fish_V3$omim_gene_perc <- go_fish_V3$omim_gene_no/length(which(universe_df[which(universe_df$omim=="Y"),"gene"]%in%goslim_anno_V3$V3))

#plotting number of genes with annotations/ fractions of genes with that annotation- ALL GENES vs omim####
allgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=gene_no), stat="identity",fill="pink") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$gene_no)+300)))+
  geom_text(aes(x = go_name,y=gene_no,label = paste0(round(gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle("Non-OMIM genes CC")
allgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=gene_no), stat="identity",fill="pink") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$gene_no)+300)))+
  geom_text(aes(x = go_name,y=gene_no,label = paste0(round(gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle("Non-OMIM MF")
allgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=gene_no), stat="identity",fill="pink") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$gene_no)+300)))+
  geom_text(aes(x = go_name,y=gene_no,label = paste0(round(gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle("Non-OMIM BP")

omimgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=omim_gene_no), stat="identity",fill="#89cff0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$omim_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omim_gene_no,label = paste0(round(omim_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle("OMIM genes CC")
omimgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=omim_gene_no), stat="identity",fill="#89cff0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$omim_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omim_gene_no,label = paste0(round(omim_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle("OMIM genes MF")
omimgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=omim_gene_no), stat="identity",fill="#89cff0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$omim_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omim_gene_no,label = paste0(round(omim_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle("OMIM genes BP")

go_fish_V3$allvsomim <- go_fish_V3$omim_gene_perc-go_fish_V3$gene_perc
go_fish_V3$allvsomimcolour <- as.character(ifelse(go_fish_V3$allvsomim>0,"#89cff0","pink"))
go_fish_V3$allvsomimcolour <- factor(go_fish_V3$allvsomimcolour ,levels=c("pink","#89cff0"))
allvsomim_CC <- ggplot(go_fish_V3[1:10,])+
  geom_bar(aes(x=go_name,y=allvsomim*100,fill=allvsomimcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("pink","#89cff0"))
allvsomim_MF <- ggplot(go_fish_V3[11:20,]) +
  geom_bar(aes(x=go_name,y=allvsomim*100,fill=allvsomimcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("pink","#89cff0"))
allvsomim_BP <- ggplot(go_fish_V3[21:30,]) +
  geom_bar(aes(x=go_name,y=allvsomim*100,fill=allvsomimcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#89cff0"))


ggarrange(
  ggarrange(allgenes_nos_BP,omimgenes_nos_BP,allvsomim_BP,ncol=3),
  ggarrange(allgenes_nos_MF,omimgenes_nos_MF,allvsomim_MF,ncol =3),
  ggarrange(allgenes_nos_CC,omimgenes_nos_CC,allvsomim_CC,ncol=3),nrow=3
)
ggsave("output/Figures/goslim_V3_allvsomim.pdf",height=30,width=30,units='cm')

#go_fish_V3 omim dom vs omim rec####
go_fish_V3$omimdom_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd"))]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$omimdom_gene_no <- lengths(go_fish_V3$omimdom_annotated_genes)
go_fish_V3$omimdom_gene_perc <- go_fish_V3$omimdom_gene_no/length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")),"gene"]%in%goslim_anno_V3$V3))

go_fish_V3$omimrec_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR"))]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$omimrec_gene_no <- lengths(go_fish_V3$omimrec_annotated_genes)
go_fish_V3$omimrec_gene_perc <- go_fish_V3$omimrec_gene_no/length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))


#plotting number of genes with annotations/ fractions of genes with that annotation- omim dom vs omim rec####
omimdomgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=omimdom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$omimdom_gene_no)+20)))+
  geom_text(aes(x = go_name,y=omimdom_gene_no,label = paste0(round(omimdom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle(paste0("Dominant OMIM genes CC   n = ",length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")),"gene"]%in%goslim_anno_V3$V3))))
omimdomgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=omimdom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$omimdom_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omimdom_gene_no,label = paste0(round(omimdom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle(paste0("Dominant OMIM genes MF   n = ",length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")),"gene"]%in%goslim_anno_V3$V3))))
omimdomgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=omimdom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$omimdom_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omimdom_gene_no,label = paste0(round(omimdom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle(paste0("Dominant OMIM genes BP   n = ",length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")),"gene"]%in%goslim_anno_V3$V3))))

omimrecgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=omimrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$omimrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omimrec_gene_no,label = paste0(round(omimrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle(paste0("Recessive OMIM genes CC   n = ", length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))))
omimrecgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=omimrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$omimrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omimrec_gene_no,label = paste0(round(omimrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle(paste0("Recessive OMIM genes MF   n = ", length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))))
omimrecgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=omimrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$omimrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=omimrec_gene_no,label = paste0(round(omimrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.03))+
  ggtitle(paste0("Recessive OMIM genes BP   n = ", length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$Inheritance_pattern=="AR"|universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))))

go_fish_V3$domomimvsrecomim <- go_fish_V3$omimdom_gene_perc-go_fish_V3$omimrec_gene_perc
go_fish_V3$domomimvsrecomimcolour <- as.character(ifelse(go_fish_V3$domomimvsrecomim>0,"dom","rec"))
go_fish_V3$domomimvsrecomimcolour <- factor(go_fish_V3$domomimvsrecomimcolour ,levels=c("dom","rec"))
domomimvsrecomim_CC <- ggplot(go_fish_V3[1:10,])+
  geom_bar(aes(x=go_name,y=domomimvsrecomim*100,fill=domomimvsrecomimcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))
domomimvsrecomim_MF <- ggplot(go_fish_V3[11:20,]) +
  geom_bar(aes(x=go_name,y=domomimvsrecomim*100,fill=domomimvsrecomimcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))
domomimvsrecomim_BP <- ggplot(go_fish_V3[21:30,]) +
  geom_bar(aes(x=go_name,y=domomimvsrecomim*100,fill=domomimvsrecomimcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))

ggarrange(
  ggarrange(omimdomgenes_nos_BP,omimrecgenes_nos_BP,domomimvsrecomim_BP,ncol=3),
  ggarrange(omimdomgenes_nos_MF,omimrecgenes_nos_MF,domomimvsrecomim_MF,ncol=3),
  ggarrange(omimdomgenes_nos_CC,omimrecgenes_nos_CC,domomimvsrecomim_CC,ncol=3),nrow=3
)
ggsave("output/Figures/goslim_V3_domomimvsrecomim.pdf",height=30,width=30,units='cm')

#go_fish_V3 lethal dom vs lethal rec####
go_fish_V3$lethaldom_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd"))]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$lethaldom_gene_no <- lengths(go_fish_V3$lethaldom_annotated_genes)
go_fish_V3$lethaldom_gene_perc <- go_fish_V3$lethaldom_gene_no/length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),"gene"]%in%goslim_anno_V3$V3))

go_fish_V3$lethalrec_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR"))]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$lethalrec_gene_no <- lengths(go_fish_V3$lethalrec_annotated_genes)
go_fish_V3$lethalrec_gene_perc <- go_fish_V3$lethalrec_gene_no/length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))


#plotting number of genes with annotations/ fractions of genes with that annotation- lethal dom vs lethal rec####
lethaldomgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=lethaldom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$lethaldom_gene_no)+5)))+
  geom_text(aes(x = go_name,y=lethaldom_gene_no,label = paste0(round(lethaldom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Dominant Lethal genes CC   n = ",length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),"gene"]%in%goslim_anno_V3$V3))))
lethaldomgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=lethaldom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$lethaldom_gene_no)+5)))+
  geom_text(aes(x = go_name,y=lethaldom_gene_no,label = paste0(round(lethaldom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Dominant Lethal genes MF   n = ",length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),"gene"]%in%goslim_anno_V3$V3))))
lethaldomgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=lethaldom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$lethaldom_gene_no)+5)))+
  geom_text(aes(x = go_name,y=lethaldom_gene_no,label = paste0(round(lethaldom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Dominant Lethal genes BP   n = ",length(which(universe_df[which(universe_df$human_lethal_B=="Y"&(universe_df$lethal_inheritance=="AD"|universe_df$lethal_inheritance=="XLd"|universe_df$lethal_inheritance=="MT,XLd")),"gene"]%in%goslim_anno_V3$V3))))

lethalrecgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=lethalrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$lethalrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=lethalrec_gene_no,label = paste0(round(lethalrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Recessive Lethal genes CC   n = ", length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))))
lethalrecgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=lethalrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$lethalrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=lethalrec_gene_no,label = paste0(round(lethalrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Recessive Lethal genes MF   n = ", length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))))
lethalrecgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=lethalrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$lethalrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=lethalrec_gene_no,label = paste0(round(lethalrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Recessive Lethal genes BP   n = ", length(which(universe_df[which(universe_df$omim=="Y"&(universe_df$lethal_inheritance=="AR"|universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="MT,AR")),"gene"]%in%goslim_anno_V3$V3))))

go_fish_V3$domlethalvsreclethal <- go_fish_V3$lethaldom_gene_perc-go_fish_V3$lethalrec_gene_perc
go_fish_V3$domlethalvsreclethalcolour <- as.character(ifelse(go_fish_V3$domlethalvsreclethal>0,"dom","rec"))
go_fish_V3$domlethalvsreclethalcolour <- factor(go_fish_V3$domlethalvsreclethalcolour ,levels=c("dom","rec"))
domlethalvsreclethal_CC <- ggplot(go_fish_V3[1:10,])+
  geom_bar(aes(x=go_name,y=domlethalvsreclethal*100,fill=domlethalvsreclethalcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))
domlethalvsreclethal_MF <- ggplot(go_fish_V3[11:20,]) +
  geom_bar(aes(x=go_name,y=domlethalvsreclethal*100,fill=domlethalvsreclethalcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))
domlethalvsreclethal_BP <- ggplot(go_fish_V3[21:30,]) +
  geom_bar(aes(x=go_name,y=domlethalvsreclethal*100,fill=domlethalvsreclethalcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))

ggarrange(
  ggarrange(lethaldomgenes_nos_BP,lethalrecgenes_nos_BP,domlethalvsreclethal_BP,ncol=3),
  ggarrange(lethaldomgenes_nos_MF,lethalrecgenes_nos_MF,domlethalvsreclethal_MF,ncol=3),
  ggarrange(lethaldomgenes_nos_CC,lethalrecgenes_nos_CC,domlethalvsreclethal_CC,ncol=3),nrow=3
)
ggsave("output/Figures/goslim_V3_domlethalvsreclethal.pdf",height=30,width=30,units='cm')









#go_fish_V3 putative lethal dom vs putative lethal rec####
go_fish_V3$putlethaldom_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y")]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$putlethaldom_gene_no <- lengths(go_fish_V3$putlethaldom_annotated_genes)
go_fish_V3$putlethaldom_gene_perc <- go_fish_V3$putlethaldom_gene_no/length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno_V3$V3))

go_fish_V3$putlethalrec_annotated_genes <- lapply(go_fish_V3$go_id,function(x) 
  unique(goslim_anno_V3$V3[which(as.character(goslim_anno_V3$V3)%in%universe_df$gene[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N")]&as.character(goslim_anno_V3$V5)==x)]))
go_fish_V3$putlethalrec_gene_no <- lengths(go_fish_V3$putlethalrec_annotated_genes)
go_fish_V3$putlethalrec_gene_perc <- go_fish_V3$putlethalrec_gene_no/length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno_V3$V3))
#plotting number of genes with annotations/ fractions of genes with that annotation- putative lethal dom vs lethal rec####
putlethaldomgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=putlethaldom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$putlethaldom_gene_no)+50)))+
  geom_text(aes(x = go_name,y=putlethaldom_gene_no,label = paste0(round(putlethaldom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Putative dominant Lethal genes CC   n = ",length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno_V3$V3))))
putlethaldomgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=putlethaldom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$putlethaldom_gene_no)+50)))+
  geom_text(aes(x = go_name,y=putlethaldom_gene_no,label = paste0(round(putlethaldom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Putative dominant Lethal genes MF   n = ",length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno_V3$V3))))
putlethaldomgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=putlethaldom_gene_no), stat="identity",fill="#008b45") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$putlethaldom_gene_no)+50)))+
  geom_text(aes(x = go_name,y=putlethaldom_gene_no,label = paste0(round(putlethaldom_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Putative dominant Lethal genes BP   n = ",length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"),"gene"]%in%goslim_anno_V3$V3))))

putlethalrecgenes_nos_CC<-ggplot(go_fish_V3[1:10,]) +
  geom_bar( aes(x=go_name, y=putlethalrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[1:10,]$putlethalrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=putlethalrec_gene_no,label = paste0(round(putlethalrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Recessive Lethal genes CC   n = ", length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno_V3$V3))))
putlethalrecgenes_nos_MF<-ggplot(go_fish_V3[11:20,]) +
  geom_bar( aes(x=go_name, y=putlethalrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[11:20,]$putlethalrec_gene_no)+50)))+
  geom_text(aes(x = go_name,y=putlethalrec_gene_no,label = paste0(round(putlethalrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Recessive Lethal genes MF   n = ", length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno_V3$V3))))
putlethalrecgenes_nos_BP<-ggplot(go_fish_V3[21:30,]) +
  geom_bar( aes(x=go_name, y=putlethalrec_gene_no), stat="identity",fill="#a020f0") +
  bar_theme()+xlab("")+
  theme(axis.text.x = element_text(angle=60, hjust=1))+scale_y_continuous(expand = c(0, 0),lim=c(0,(max(go_fish_V3[21:30,]$putlethalrec_gene_no)+100)))+
  geom_text(aes(x = go_name,y=putlethalrec_gene_no,label = paste0(round(putlethalrec_gene_perc*100,0),"%")),size=3,position = position_stack(vjust = 1.1))+
  ggtitle(paste0("Recessive Lethal genes BP   n = ",length(which(universe_df[which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="N"),"gene"]%in%goslim_anno_V3$V3))))

go_fish_V3$putdomlethalvsreclethal <- go_fish_V3$putlethaldom_gene_perc-go_fish_V3$putlethalrec_gene_perc
go_fish_V3$putdomlethalvsreclethalcolour <- as.character(ifelse(go_fish_V3$putdomlethalvsreclethal>0,"dom","rec"))
go_fish_V3$putdomlethalvsreclethalcolour <- factor(go_fish_V3$putdomlethalvsreclethalcolour ,levels=c("dom","rec"))
putdomlethalvsreclethal_CC <- ggplot(go_fish_V3[1:10,])+
  geom_bar(aes(x=go_name,y=putdomlethalvsreclethal*100,fill=putdomlethalvsreclethalcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))
putdomlethalvsreclethal_MF <- ggplot(go_fish_V3[11:20,]) +
  geom_bar(aes(x=go_name,y=putdomlethalvsreclethal*100,fill=putdomlethalvsreclethalcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))
putdomlethalvsreclethal_BP <- ggplot(go_fish_V3[21:30,]) +
  geom_bar(aes(x=go_name,y=putdomlethalvsreclethal*100,fill=putdomlethalvsreclethalcolour),stat="identity")+bar_theme()+geom_hline(yintercept=0)+
  theme(axis.text.x = element_text(angle=60, hjust=1),axis.line.x=element_blank(),legend.position = "none")+scale_fill_manual(values=c("#008b45","#a020f0"))

ggarrange(
  ggarrange(putlethaldomgenes_nos_BP,putlethalrecgenes_nos_BP,putdomlethalvsreclethal_BP,ncol=3),
  ggarrange(putlethaldomgenes_nos_MF,putlethalrecgenes_nos_MF,putdomlethalvsreclethal_MF,ncol=3),
  ggarrange(putlethaldomgenes_nos_CC,putlethalrecgenes_nos_CC,putdomlethalvsreclethal_CC,ncol=3),nrow=3
)
ggsave("output/Figures/goslim_V3_putdomlethalvsreclethal.pdf",height=30,width=30,units='cm')











