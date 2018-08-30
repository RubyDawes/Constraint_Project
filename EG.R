############yeast##########
yeast_viable <-read.table("Gene_lists/EG/yeast/viable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable <-read.table("Gene_lists/EG/yeast/inviable_annotations.txt",quote="",header=FALSE,comment.char="!",fill=TRUE,sep="\t")
yeast_inviable<-yeast_inviable[which(yeast_inviable$V6=="null"),]
yeast_viable<-yeast_viable[which(yeast_viable$V6=="null"),]
length(unique(yeast_inviable$V1))+length(unique(yeast_viable$V1))
############pie of yeast lethal genes##############
yeast_lethal_pie <- data.frame(group=c("lethal", "non-lethal"),
                             value=c(
                               length(unique(yeast_inviable$V1)),
                               length(unique(yeast_viable$V1))))
yeast_lethal_pie$group <- factor(yeast_lethal_pie$group, levels = yeast_lethal_pie$group)
c <- ggplot(yeast_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("black","steelblue3"))
c<- c+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/yeastpie.pdf",height=18, width=18, units='cm')

############fly##########
#flyerror<-read.table("Gene_lists/EG/fly/allele_phenotypic_data_fb_2018_03.tsv",header=FALSE,comment.char="#",fill=TRUE,sep="\t")

fly<-read.table("Gene_lists/EG/fly/allele_phenotypic_data_fb_2018_03.tsv",quote="",header=FALSE,comment.char="#",fill=TRUE,sep="\t")
fly <- fly[which(!grepl('\\\\',fly$V5)&!grepl('\\::',fly$V5)),]
#require(dplyr) 
#weird1<-anti_join(fly,flyerror)
#weird2<-anti_join(flyerror,fly)
#length(unique(weird1$V5[which(grepl("lethal",weird1$V3))]))
#length(unique(weird2$V5[which(grepl("lethal",weird2$V3))]))
#a<-which(grepl("lethal",weird1$V3))
#b<-which(grepl("lethal",weird2$V3))
#rlyweird<-weird1[a[which(!a%in%b)],]
#length(unique(fly$V5))
flygenes<-read.table("Gene_lists/EG/fly/fbgn_NAseq_Uniprot_fb_2018_04.tsv",quote="",header=FALSE,comment.char="#",fill=TRUE,sep="\t")

length(unique(fly$V5))
fly_lethal <-fly[which(grepl("lethal",fly$V3)&grepl("stage",fly$V3)),]
length(unique(fly_lethal$V5))

#fly_lof <- read.table("Gene_lists/EG/fly/FlyBase_IDs_lof.txt",header=FALSE,comment.char="#",fill=TRUE,sep="\t")
#fly_lof <- fly[which(fly$V2%in%fly_lof$V1),]
#length(unique(fly_lof$V5))

#fly_lof_lethal <- fly_lof[which(grepl("lethal",fly_lof$V3)),]
#length(unique(fly_lof_lethal$V5))

############pie of fly lethal genes##############
fly_lethal_pie <- data.frame(group=c("lethal", "non-lethal"),
                               value=c(
                                 length(unique(fly_lethal$V5)),
                                 length(which(!unique(fly$V5)%in%unique(fly_lethal$V5)))))
fly_lethal_pie$group <- factor(fly_lethal_pie$group, levels = fly_lethal_pie$group)
c <- ggplot(fly_lethal_pie, aes(x="", y=value, fill=group))
c<- c+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
c<- c+scale_fill_manual(values=c("black","steelblue3"))
c<- c+geom_text(aes(y = value,label = percent(value/sum(value))),size=6,position = position_stack(vjust = 0.5))
ggsave("Analysis/Sandra_Figures/Figs/flypie.pdf",height=18, width=18, units='cm')

############fish##########
fish<-read.table("Gene_lists/EG/FISH/phenoGeneCleanData_fish_2018.08.26.txt",comment.char="#",fill=TRUE,sep="\t")
fish_info<-read.table("Gene_lists/EG/FISH/features-affected-genes_2018.08.26.txt",comment.char="#",fill=TRUE,sep="\t")



length(unique(fish$V2))
fish_lethal<-fish[which(grepl("lethal",fish$V11)),]
length(unique(fish_lethal$V2))

length(which(grepl("whole organism",fish$V10)))