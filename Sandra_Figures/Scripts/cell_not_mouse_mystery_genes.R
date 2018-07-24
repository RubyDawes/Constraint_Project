cell_not_mouse <- universe_df[which(universe_df$cell_essential=="Y"&universe_df$lethal_mouse=="N"),]

length(which(cell_not_mouse$cell_essential_hits=="3"))

high_phens <- c("adipose tissue phenotype","behavior/neurological phenotype","cardiovascular system phenotype","cellular phenotype",
                "craniofacial phenotype","digestive/alimentary phenotype","embryo phenotype","endocrine/exocrine gland phenotype",
                "growth/size/body region phenotype","hearing/vestibular/ear phenotype","hematopoietic system phenotype","homeostasis/metabolism phenotype",
                "immune system phenotype","integument phenotype","limbs/digits/tail phenotype","liver/biliary system phenotype","mortality/aging",
                "muscle phenotype","neoplasm","nervous system phenotype","normal phenotype","pigmentation phenotype","renal/urinary system phenotype",
                "reproductive system phenotype","respiratory system phenotype","skeleton phenotype","taste/olfaction phenotype","vision/eye phenotype"
)

for (i in 1:length(high_phens)) {
  cell_not_mouse[,53+i] <- ifelse(grepl(high_phens[i],cell_not_mouse$high_MP_phen),1,0)
}
rm(i)
colnames(cell_not_mouse)[54:ncol(cell_not_mouse)] <- high_phens

high_phens <- c("adipose tissue ","behavior/neurological ","cardiovascular system ","cellular ",
                "craniofacial ","digestive/alimentary ","embryo ","endocrine/exocrine gland ",
                "growth/size/body region ","hearing/vestibular/ear ","hematopoietic system ","homeostasis/metabolism ",
                "immune system ","integument ","limbs/digits/tail ","liver/biliary system ","mortality/aging",
                "muscle ","neoplasm","nervous system ","normal ","pigmentation ","renal/urinary system ",
                "reproductive system ","respiratory system ","skeleton ","taste/olfaction ","vision/eye "
)
#for each high phen how many cell_not_mouse genes had that phen category
percs_phen <- c()
for (i in 1:length(high_phens)){
  percs_phen[i] <- length(which(cell_not_mouse[53+i]==1))/length(which(!is.na(cell_not_mouse[53+i])))
}
rm(i)
a<- data.frame(high_phens,percs_phen)
b <- ggplot(melt(a),aes(x=reorder(high_phens,value),value))+geom_bar(aes(fill=variable),position="dodge",stat="identity")
b<- b+coord_flip()+bar_twosided_theme()+geom_hline(yintercept=0)+labs(y="",x="High-level phenotype")
b <- b+scale_y_continuous(breaks = pretty(a$perc_diff, n = 6))+ggtitle("Mouse phenotypes of Cell essential mouse non-lethal genes")
png("output/figures/cellessential_phens.png",width=1800,height=1000,type="quartz",res=150,bg = "transparent")
b
dev.off()


