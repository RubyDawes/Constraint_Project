#DNA/RNA synthesis repair replication
cat1<-universe_df$gene[which(grepl("DNA repair",universe_df$go_names)|
                         grepl("DNA replication",universe_df$go_names)|
                         grepl("DNA biosynthetic process",universe_df$go_names)|
                         grepl("chromosome segregation",universe_df$go_names)|
                         grepl("RNA biosynthetic process",universe_df$go_names)|
                         grepl("RNA repair",universe_df$go_names)|
                         grepl("RNA replication",universe_df$go_names)|
                         grepl("RNA processing",universe_df$go_names)|
                         grepl("RNA catabolic process",universe_df$go_names)|
                         grepl("DNA catabolic process",universe_df$go_names)|
                         grepl("RNA polymerase",universe_df$go_names)|
                         grepl("RNA splicing",universe_df$go_names)|
                         grepl("RNA metabolic process",universe_df$go_names)|
                         grepl("ribonuclease",universe_df$go_names)|
                         grepl("ribonucleoside",universe_df$go_names)|
                         grepl("chromatin organization",universe_df$go_names)|
                         grepl("spindle assembly",universe_df$go_names)|
                         grepl("strand break repair",universe_df$go_names)|
                         grepl("transcription",universe_df$go_names)
                       )]

#Protein synthesis and PTM
cat2<-universe_df$gene[which(grepl("regulation of protein folding",universe_df$go_names)|
                         grepl("protein folding",universe_df$go_names)|
                         grepl("post-translational protein modification",universe_df$go_names)|
                           grepl("phosphorylation",universe_df$go_names)|
                           grepl("methyltransferase",universe_df$go_names)|
                         grepl("translation",universe_df$go_names)|
                         grepl("tRNA",universe_df$go_names)
                      )]


#Protein degradation/ Proteases
cat3<-universe_df$gene[which(grepl("peptidase activity",universe_df$go_names)|
                               grepl("ubiquitin",universe_df$go_names)|
                         grepl("protein catabolic process",universe_df$go_names)
)]

#Lipid biosynthesis/energy metabolism/redox
cat4<-universe_df$gene[which(grepl("lipid biosynthetic process",universe_df$go_names)|
                         grepl("generation of precursor metabolites and energy",universe_df$go_names)|
                         grepl("oxidoreductase activity",universe_df$go_names)|
                         grepl("cell redox homeostasis",universe_df$go_names)|
                         grepl("regulation of cellular metabolic process",universe_df$go_names)|
                         grepl("lipid metabolic process",universe_df$go_names)|
                         grepl("lipid catabolic process",universe_df$go_names)|
                           grepl("oligosaccharide",universe_df$go_names)|
                           grepl("mannose",universe_df$go_names)|
                           grepl("glycosylation",universe_df$go_names)|
                           grepl("glycan",universe_df$go_names)|
                         grepl("fatty acid",universe_df$go_names)
)]

#Mitochondrial
cat5<-universe_df$gene[which(grepl("mitochondria",universe_df$go_names)
)]

  
#Organelle biogenesis/membrane traffic
cat6<-universe_df$gene[which(grepl("regulation of cellular component biogenesis",universe_df$go_names)|
                         grepl("cellular component organization or biogenesis",universe_df$go_names)|
                         grepl("maintenance of location in cell",universe_df$go_names)|
                         grepl("organelle organization",universe_df$go_names)|
                         grepl("organelle assembly",universe_df$go_names)|
                         grepl("transport vesicle",universe_df$go_names)|
                         grepl("GTPase",universe_df$go_names)|
                         grepl("golgi organisation",universe_df$go_names)|
                         grepl("vesicle",universe_df$go_names)|
                         grepl("transport",universe_df$go_names))]


#Transporters+channels
cat7<-universe_df$gene[which(grepl("transporter",universe_df$go_names)|
                         grepl("channel",universe_df$go_names)|
                         grepl("transmembrane movement",universe_df$go_names)
                         )]


#Adhesion/ base membrane/ ECM
cat8<-universe_df$gene[which(grepl("cell adhesion",universe_df$go_names)|
                         grepl("regulation of cell adhesion",universe_df$go_names)|
                         grepl("basement membrane",universe_df$go_names)|
                         grepl("basement membrane organization",universe_df$go_names)|
                         grepl("basement membrane assembly",universe_df$go_names)|
                         grepl("extracellular matrix organization",universe_df$go_names)|
                         grepl("extracellular matrix",universe_df$go_names)
)]

#Cytoskeletal proteins + motors
cat9<-universe_df$gene[which(grepl("cytoskeleton",universe_df$go_names)|
                         grepl("regulation of cellular component movement",universe_df$go_names)|
                         grepl("actin filament-based process",universe_df$go_names)|
                         grepl("intermediate filament-based process",universe_df$go_names)|
                         grepl("microtubule-based process",universe_df$go_names)|
                         grepl("motor activity",universe_df$go_names)
)]

#Signalling and development
cat10<-universe_df$gene[which(grepl("signaling",universe_df$go_names)|
                              grepl("signal transduction",universe_df$go_names)|
                              grepl("receptor",universe_df$go_names)|
                                grepl("regulation of gene expression",universe_df$go_names)|
                              grepl("development",universe_df$go_names))]

#Immune system 
cat11<-universe_df$gene[which(grepl("immune",universe_df$go_names)|
                                grepl("cytokine",universe_df$go_names)|
                                grepl("development",universe_df$go_names))]

#cell cycle
cat12<-universe_df$gene[which(grepl("cell cycle",universe_df$go_names))]



unsorted<- universe_df[which(lengths(universe_df$go_names)>0&!universe_df$gene%in%union(cat1,union(cat2,union(cat3,union(cat4,union(cat5,union(cat6,union(cat7,union(cat8,union(cat9,union(cat10,union(cat11,cat12)))))))))))),c(2,33)]


#annotation/ sortng rate
all_annotated<-union(cat1,union(cat2,union(cat3,union(cat4,union(cat5,union(cat6,union(cat7,union(cat8,union(cat9,union(cat10,union(cat11,cat12)))))))))))
sorting_rate <- data.frame(group=c("unannotated","annotated,unsorted","sorted"),
                              value=c(length(which(lengths(universe_df$go_names)==0)),
                                       length(which(lengths(universe_df$go_names)>0&!universe_df$gene%in%all_annotated)),
                                       length(all_annotated))
                              )
sorting_rate$group <- factor(sorting_rate$group, levels =sorting_rate$group)

a <- ggplot(sorting_rate, aes(x="", y=value, fill=group))
a<- a+geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)
a<- a+scale_fill_manual(values=c("slategrey","snow4","springgreen3"))#+theme(legend.position="none")
a<- a+geom_text(aes(y = value,label = value),size=6,position = position_stack(vjust = 0.5))
a<- a+ggtitle(paste("Sorting of all protein coding genes into custom functional categories"))+theme(plot.title = element_text(size = 12, face = "bold"),legend.position="none",legend.text=element_text(size=12))

#numbers in each category
length(which(universe_df$cell_essential=="Y"&universe_df$cat1=="Y"))/length(which(universe_df$cat1=="Y"&!is.na(universe_df$cell_essential)))
length(which(universe_df$lethal_mouse=="Y"&universe_df$cat1=="Y"))/length(which(universe_df$cat1=="Y"&!is.na(universe_df$lethal_mouse)))


#annotating universe_df with custom functional categories
universe_df$fun_cat <- rep(NA,length(universe_df$gene))
universe_df$cat1 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat1,"Y","N")))
universe_df$cat2 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat2,"Y","N")))
universe_df$cat3 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat3,"Y","N")))
universe_df$cat4 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat4,"Y","N")))
universe_df$cat5 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat5,"Y","N")))
universe_df$cat6 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat6,"Y","N")))
universe_df$cat7 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat7,"Y","N")))
universe_df$cat8 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat8,"Y","N")))
universe_df$cat9 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat9,"Y","N")))
universe_df$cat10 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat10,"Y","N")))
universe_df$cat11 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat11,"Y","N")))
universe_df$cat12 <- unlist(lapply(universe_df$gene, function(x) ifelse(x%in%cat12,"Y","N")))

universe_df$fun_cat[which(universe_df$cat1=="Y")]<-"DNA/RNA synthesis, repair and replication"
universe_df$fun_cat[which(universe_df$cat2 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat2 == "Y")]),"Protein synthesis and PTM",paste(universe_df$fun_cat[which(universe_df$cat2 == "Y")],"Protein synthesis and PTM",sep="; "))
universe_df$fun_cat[which(universe_df$cat3 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat3 == "Y")]),"Protein degradation/ Proteases",paste(universe_df$fun_cat[which(universe_df$cat3 == "Y")],"Protein degradation/ Proteases",sep="; "))
universe_df$fun_cat[which(universe_df$cat4 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat4 == "Y")]),"Lipid biosynthesis/energy metabolism/redox",paste(universe_df$fun_cat[which(universe_df$cat4 == "Y")],"Lipid biosynthesis/energy metabolism/redox",sep="; "))
universe_df$fun_cat[which(universe_df$cat5 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat5 == "Y")]),"Mitochondrial",paste(universe_df$fun_cat[which(universe_df$cat5 == "Y")],"Mitochondrial",sep="; "))
universe_df$fun_cat[which(universe_df$cat6 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat6 == "Y")]),"Organelle biogenesis/membrane traffic",paste(universe_df$fun_cat[which(universe_df$cat6 == "Y")],"Organelle biogenesis/membrane traffic",sep="; "))
universe_df$fun_cat[which(universe_df$cat7 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat7 == "Y")]),"Transporters+channels",paste(universe_df$fun_cat[which(universe_df$cat7 == "Y")],"Transporters+channels",sep="; "))
universe_df$fun_cat[which(universe_df$cat8 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat8 == "Y")]),"Adhesion/ base membrane/ ECM",paste(universe_df$fun_cat[which(universe_df$cat8 == "Y")],"Adhesion/ base membrane/ ECM",sep="; "))
universe_df$fun_cat[which(universe_df$cat9 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat9 == "Y")]),"Cytoskeletal proteins + motors",paste(universe_df$fun_cat[which(universe_df$cat9 == "Y")],"Cytoskeletal proteins + motors",sep="; "))
universe_df$fun_cat[which(universe_df$cat10 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat10 == "Y")]),"Signalling and development",paste(universe_df$fun_cat[which(universe_df$cat10 == "Y")],"Signalling and development",sep="; "))
universe_df$fun_cat[which(universe_df$cat11 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat11 == "Y")]),"Immune system ",paste(universe_df$fun_cat[which(universe_df$cat11 == "Y")],"Immune system ",sep="; "))
universe_df$fun_cat[which(universe_df$cat12 == "Y")]<- ifelse(is.na(universe_df$fun_cat[which(universe_df$cat12 == "Y")]),"cell cycle",paste(universe_df$fun_cat[which(universe_df$cat12 == "Y")],"cell cycle",sep="; "))

universe_df<-universe_df[,-c(42:53)]

