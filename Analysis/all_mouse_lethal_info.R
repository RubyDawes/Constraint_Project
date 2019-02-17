#preparing df of all mouse info from MGI####
mgi_df <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenoGenoMP.xlsx")

mgi_df <- mgi_df[,c(1,2,6,4)]

#removing any rows with a comma in the MGI ID column - indicates a multi-gene knockout (we are only interested in single gene knockouts)
mgi_df<- mgi_df[which(!grepl(pattern = "\\,", x = mgi_df$MGI_ID, ignore.case = FALSE)),]

#zygosity column
mgi_df$zygosity <- ifelse(grepl("\\|",mgi_df$`Allele.Symbol(s)`),"heterozygous","homozygous")
  
#getting allele info for heterozygous mice
#extracting only non-wildtype allele
mgi_df$`Allele.Symbol(s)`[which(mgi_df$zygosity=="heterozygous")] <- strsplit(mgi_df$`Allele.Symbol(s)`[which(mgi_df$zygosity=="heterozygous")],"\\|")
mgi_df$`Allele.Symbol(s)`[which(mgi_df$zygosity=="heterozygous")]  <- lapply(which(mgi_df$zygosity=="heterozygous"), function(x) mgi_df$`Allele.Symbol(s)`[[x]][which(grepl("\\+",mgi_df$`Allele.Symbol(s)`[[x]])==FALSE)])

mgi_df$zygosity[which(lengths(mgi_df$`Allele.Symbol(s)`)>1)]<- "compound heterozygous"

mgi_lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
mgi_df$is_lethal <- ifelse(mgi_df$MP_ID%in%mgi_lethalphens$MP.id,"Y","N")


#reading in allele data 
allele <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenotypicAllele.xlsx")

mgi_df$allele_info <- lapply(mgi_df$`Allele.Symbol(s)`,function(x) vlookup(x,allele,result_column="Allele.Attribute",lookup_column="Allele.Symbol"))
mgi_df$allele_type <- lapply(mgi_df$`Allele.Symbol(s)`,function(x) vlookup(x,allele,result_column="Allele.Type",lookup_column="Allele.Symbol"))
mgi_df$allele_ID <- lapply(mgi_df$`Allele.Symbol(s)`,function(x) vlookup(x,allele,result_column="MGI.Allele.Accession.ID",lookup_column="Allele.Symbol"))

# add on gene names and human ortholog names
mgi_names <- read.xlsx("Gene_lists/MGI_MPs/HMD_HumanPhenotype.xlsx")
mgi_names$MGI.Marker.Accession.ID <- gsub(" ", "", mgi_names$MGI.Marker.Accession.ID, fixed = TRUE)
mgi_df$mouse_symbol <- vlookup(mgi_df$MGI_ID,mgi_names,result_column="Mouse.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi_df$human_symbol <- vlookup(mgi_df$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi_df$human_symbol <- checkGeneSymbols(mgi_df$human_symbol,unmapped.as.na=FALSE)[[3]]

mgi_df$ko <- lapply(mgi_df$allele_info,function(x)
  ifelse(grepl("Null/knockout|Hypomorph",x),"Y","N"))
mgi_df$ko <- lapply(mgi_df$ko,unique)

write.xlsx(mgi_df,"output/spreadsheets/MGI_all_info.xlsx")
save(mgi_df, file="output/Data/MGI_all_data.rda", compress="bzip2")

#analysis ####
load("output/Data/MGI_all_data.rda")

#converting into per-gene info

mgi_bd <- data.frame(gene = unique(mgi_df$human_symbol[which(!is.na(mgi_df$human_symbol))]))
mgi_bd$MGI_ID = lapply(mgi_bd$gene, function(x) unique(mgi_df$MGI_ID[which(mgi_df$human_symbol==x)]))

mgi_homko <- mgi_df[which(mgi_df$zygosity=="homozygous"&mgi_df$ko=="Y"),]
mgi_bd$hom_ko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_homko$MGI_ID,"Y","N"))
mgi_bd$hom_ko <- unlist(lapply(mgi_bd$hom_ko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$hom_ko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_homko$is_lethal[which(mgi_homko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

mgi_hetko <- mgi_df[which(mgi_df$zygosity=="heterozygous"&mgi_df$ko=="Y"),]
mgi_bd$het_ko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_hetko$MGI_ID,"Y","N"))
mgi_bd$het_ko <- unlist(lapply(mgi_bd$het_ko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$het_ko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_hetko$is_lethal[which(mgi_hetko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

mgi_comphet2ko <- mgi_df[which(mgi_df$zygosity=="compound heterozygous"&mgi_df$ko=="Y"),]
mgi_bd$comphet_2ko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_comphet2ko$MGI_ID,"Y","N"))
mgi_bd$comphet_2ko <- unlist(lapply(mgi_bd$comphet_2ko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$comphet_2ko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_comphet2ko$is_lethal[which(mgi_comphet2ko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

mgi_comphet1ko <- mgi_df[which(mgi_df$zygosity=="compound heterozygous"&lengths(mgi_df$ko)>1),]
mgi_bd$comphet_1ko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_comphet1ko$MGI_ID,"Y","N"))
mgi_bd$comphet_1ko <- unlist(lapply(mgi_bd$comphet_1ko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$comphet_1ko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_comphet1ko$is_lethal[which(mgi_comphet1ko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

mgi_comphet_nonko <- mgi_df[which(mgi_df$zygosity=="compound heterozygous"&mgi_df$ko=="N"),]
mgi_bd$comphet_nonko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_comphet_nonko$MGI_ID,"Y","N"))
mgi_bd$comphet_nonko <- unlist(lapply(mgi_bd$comphet_nonko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$comphet_nonko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_comphet_nonko$is_lethal[which(mgi_comphet_nonko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

mgi_hetnonko <- mgi_df[which(mgi_df$zygosity=="heterozygous"&mgi_df$ko=="N"),]
mgi_bd$het_nonko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_hetnonko$MGI_ID,"Y","N"))
mgi_bd$het_nonko <- unlist(lapply(mgi_bd$het_nonko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$het_nonko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_hetnonko$is_lethal[which(mgi_hetnonko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

mgi_homnonko <- mgi_df[which(mgi_df$zygosity=="homozygous"&mgi_df$ko=="N"),]
mgi_bd$hom_nonko <- lapply(mgi_bd$MGI_ID,function(x) ifelse(x%in%mgi_homnonko$MGI_ID,"Y","N"))
mgi_bd$hom_nonko <- unlist(lapply(mgi_bd$hom_nonko, function(x) ifelse(length(unique(x))>1,"Y",unique(x))))
mgi_bd$hom_nonko_lethal <- unlist(lapply(1:length(mgi_bd$gene),function(x) ifelse("Y"%in%mgi_homnonko$is_lethal[which(mgi_homnonko$MGI_ID%in%mgi_bd$MGI_ID[[x]])],"Y","N")))

#summary for each gene
mgi_bd$lethality <- rep(NA,length(mgi_bd$gene))
mgi_bd$lethality[which(mgi_bd$hom_ko_lethal=="Y")] <- "Hom/z KO lethal"
mgi_bd$lethality[which(mgi_bd$het_ko_lethal=="Y")] <- ifelse(is.na(mgi_bd$lethality[which(mgi_bd$het_ko_lethal=="Y")]),
                                                             " Het/z KO lethal",paste(mgi_bd$lethality[which(mgi_bd$het_ko_lethal=="Y")]," Het/z KO lethal",sep=","))
mgi_bd$lethality[which(mgi_bd$comphet_2ko_lethal=="Y")] <- ifelse(is.na(mgi_bd$lethality[which(mgi_bd$comphet_2ko_lethal=="Y")]),
                                                             " Comp Het/z  2 KO lethal",paste(mgi_bd$lethality[which(mgi_bd$comphet_2ko_lethal=="Y")]," Comp Het/z  2 KO lethal",sep=","))
mgi_bd$lethality[which(mgi_bd$comphet_1ko_lethal=="Y")] <- ifelse(is.na(mgi_bd$lethality[which(mgi_bd$comphet_1ko_lethal=="Y")]),
                                                                  " Comp Het/z  mixed alleles lethal",paste(mgi_bd$lethality[which(mgi_bd$comphet_1ko_lethal=="Y")]," Comp Het/z  mixed alleles lethal",sep=","))
mgi_bd$lethality[which(mgi_bd$comphet_nonko_lethal=="Y")] <- ifelse(is.na(mgi_bd$lethality[which(mgi_bd$comphet_nonko_lethal=="Y")]),
                                                                  " Comp Het/z  non-KO lethal",paste(mgi_bd$lethality[which(mgi_bd$comphet_nonko_lethal=="Y")]," Comp Het/z  non-KO lethal",sep=","))
mgi_bd$lethality[which(mgi_bd$het_nonko_lethal=="Y")] <- ifelse(is.na(mgi_bd$lethality[which(mgi_bd$het_nonko_lethal=="Y")]),
                                                                    " Het/z  non-KO lethal",paste(mgi_bd$lethality[which(mgi_bd$het_nonko_lethal=="Y")]," Het/z  non-KO lethal",sep=","))
mgi_bd$lethality[which(mgi_bd$hom_nonko_lethal=="Y")] <- ifelse(is.na(mgi_bd$lethality[which(mgi_bd$hom_nonko_lethal=="Y")]),
                                                                " Hom/z non-KO lethal",paste(mgi_bd$lethality[which(mgi_bd$hom_nonko_lethal=="Y")]," Hom/z  non-KO lethal",sep=","))

write.xlsx(mgi_bd,"output/spreadsheets/MGI_breakdown.xlsx")
save(mgi_bd, file="output/Data/MGI_breakdown.rda", compress="bzip2")

load("output/Data/MGI_breakdown.rda")

weirdlethal <- mgi_bd$gene[which(!grepl("Hom/z KO lethal",mgi_bd$lethality)&!is.na(mgi_bd$lethality))]
#432 mgi lethal non hom/z ko genes- how many IMPC lethal? 
weirdlethal_nonimpc <- weirdlethal[which(!weirdlethal%in%universe_df$gene[which(universe_df$lethal_IMPC=="Y")])]


mgi_df$homz_ko <- ifelse((mgi_df$zygosity=="homozygous"&mgi_df$ko=="Y"),"Y","N")
#making spreadsheet with info on 385 weird lethal non-IMPC genes
mgi_wl <- data.frame(gene = weirdlethal_nonimpc)
mgi_wl$MGI_ID <- mgi_bd$MGI_ID[which(mgi_bd$gene%in%mgi_wl$gene)]
mgi_wl$lethal_mgi_info <- vlookup(mgi_wl$gene,mgi_bd, result_column = "lethality",lookup_column ="gene" )

mgi_wl$allelic_composition <- lapply(mgi_wl$gene,function(x) unlist(mgi_df$Allelic.Composition[which(mgi_df$human_symbol == x & mgi_df$is_lethal=="Y" & mgi_df$homz_ko=="N")]))
mgi_wl$allele_ID <- lapply(mgi_wl$gene,function(x) unlist(mgi_df$allele_ID[which(mgi_df$human_symbol == x & mgi_df$is_lethal=="Y" & mgi_df$homz_ko=="N")]))
mgi_wl$allele_info <- lapply(mgi_wl$gene, function(x) paste(mgi_df$allele_info[which(mgi_df$human_symbol == x & mgi_df$is_lethal=="Y" & mgi_df$homz_ko=="N")],
                                                            mgi_df$allele_type[which(mgi_df$human_symbol == x & mgi_df$is_lethal=="Y" & mgi_df$homz_ko=="N")],sep="|"))

mgi_wl$MP_ID <- lapply(mgi_wl$gene,function(x) unlist(mgi_df$MP_ID[which(mgi_df$human_symbol == x & mgi_df$is_lethal=="Y" & mgi_df$homz_ko=="N")]))

save(mgi_wl, file = "output/Data/MGI_lethalgenes_non_recessivenull.rda", compress="bzip2")
write.xlsx(mgi_wl, "output/spreadsheets/MGI_lethal_nothomozygousko_genes.xlsx")


#pie of 385 weird lethal genes
pie <- data.frame(group=c("non-OMIM", "OMIM","Human Lethal"),
                               value=c(
                                 length(which(universe_df$gene%in%weirdlethal_nonimpc&is.na(universe_df$omim))),
                                 length(which(universe_df$gene%in%weirdlethal_nonimpc&universe_df$omim=="Y"&universe_df$human_lethal_B=="N")),
                                 length(which(universe_df$gene%in%weirdlethal_nonimpc&universe_df$human_lethal_B=="Y"))
                                 ))
pie$group <- factor(pie$group, levels = pie$group)
a <- ggplot(pie, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+blank_theme()+coord_polar(theta="y",direction=-1)+
  scale_fill_manual(values=c("steelblue3","black","#8B1A1A"))+theme(plot.title=element_text(size=10,face="plain"))+
  geom_text(aes(y = value,label = paste0(value,"\n",percent(value/sum(value))),colour=group),size=3,position = position_stack(vjust = 0.5),fontface="bold")+
  scale_colour_manual(values=c("black","white","white"))+ggtitle(paste0("Protein coding genes with \n lethal, non-homozygous \n KO phenotype in MGI \n n = ",sum(pie[,2])))
ggsave("pie_380extralethal.pdf",height=8, width=14, units='cm')


