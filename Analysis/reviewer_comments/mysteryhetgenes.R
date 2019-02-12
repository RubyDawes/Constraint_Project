# do het/z lethal mouse genes hhave pLI = 1?####
het_lethal <- universe_df[which(universe_df$lethal_het_mouse=="Y"),]

write.xlsx(het_lethal,"het_lethal.xlsx",append=TRUE)

#getting homozygous and hetereozygous alleles + phenotypes from MGI####
mgi_df <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenoGenoMP.xlsx")
mgi_df <- mgi_df[,c(1,2,6,4)]
mgi_df<- mgi_df[which(!grepl(pattern = "\\,", x = mgi_df$MGI_ID, ignore.case = FALSE)),]


#extracting heterozygous alleles
# keeping only rows with data on heterozygous KO mice (by looking for a "|" in the allele symbol column)
het <- mgi_df[which(grepl("\\+",mgi_df$`Allele.Symbol(s)`)==TRUE),]

#extracting only non-wildtype allele
het$allele <- strsplit(het$`Allele.Symbol(s)`,"\\|")
het$allele <- lapply(1:length(het$allele), function(x) het$allele[[x]][which(grepl("\\+",het$allele[[x]])==FALSE)])
het<- het[-which(lengths(het$allele)!=1),]
het$allele <- unlist(het$allele)

#reading in allele data so we can filter to only include null/knockout and hypomorph alleles. 
allele <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenotypicAllele.xlsx")
het$allele_info <- vlookup(het$allele,allele,result_column="Allele.Attribute",lookup_column="Allele.Symbol")
het <- het[which(grepl(pattern="Null/knockout",x=het$allele_info,ignore.case=FALSE)|grepl(pattern="Hypomorph",x=het$allele_info,ignore.case=TRUE)),]

het$allele_type <- vlookup(het$allele,allele,result_column="Allele.Type",lookup_column="Allele.Symbol")
het$allele_ID <- vlookup(het$allele,allele,result_column="MGI.Allele.Accession.ID",lookup_column="Allele.Symbol")
het$allele_pmID <- vlookup(het$allele,allele,result_column="PubMed.ID.for.original.reference",lookup_column="Allele.Symbol")

# add on gene names and human ortholog names
mgi_names <- read.xlsx("Gene_lists/MGI_MPs/HMD_HumanPhenotype.xlsx")
mgi_names$MGI.Marker.Accession.ID <- gsub(" ", "", mgi_names$MGI.Marker.Accession.ID, fixed = TRUE)
het$mouse_symbol <- vlookup(het$MGI_ID,mgi_names,result_column="Mouse.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
het$human_symbol <- vlookup(het$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
het$human_symbol <- checkGeneSymbols(het$human_symbol,unmapped.as.na=FALSE)[[3]]


#extracting homozygous alleles from mgi_df
mgi_df$allele_info <- vlookup(mgi_df$`Allele.Symbol(s)`,allele,result_column="Allele.Attribute",lookup_column="Allele.Symbol")
mgi_df <- mgi_df[which(grepl(pattern="Null/knockout",x=mgi_df$allele_info,ignore.case=FALSE)|grepl(pattern="Hypomorph",x=mgi_df$allele_info,ignore.case=TRUE)),]

mgi_df$allele_type <- vlookup(mgi_df$`Allele.Symbol(s)`,allele,result_column="Allele.Type",lookup_column="Allele.Symbol")
mgi_df$allele_ID <- vlookup(mgi_df$`Allele.Symbol(s)`,allele,result_column="MGI.Allele.Accession.ID",lookup_column="Allele.Symbol")
mgi_df$allele_pmID <- vlookup(mgi_df$`Allele.Symbol(s)`,allele,result_column="PubMed.ID.for.original.reference",lookup_column="Allele.Symbol")

# add on gene names and human ortholog names
mgi_df$mouse_symbol <- vlookup(mgi_df$MGI_ID,mgi_names,result_column="Mouse.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi_df$human_symbol <- vlookup(mgi_df$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi_df$human_symbol <- checkGeneSymbols(mgi_df$human_symbol,unmapped.as.na=FALSE)[[3]]

#adding on MP term definitions
mgi_allphens <- read.xlsx("Gene_lists/MGI_MPs/VOC_MammalianPhenotype.xlsx",startRow=1,colNames=FALSE)
names(mgi_allphens)<- c("MP.id","phenotype","definition")
mgi_df$MP_DEF <- vlookup(mgi_df$MP_ID, mgi_allphens,result_column = "phenotype",lookup_column = "MP.id")
het$MP_DEF <- vlookup(het$MP_ID, mgi_allphens,result_column = "phenotype",lookup_column = "MP.id")






#putting into database of variants and phenotypes homozygous and heterozygous####
naughties <- universe_df[which(universe_df$lethal_het_mouse=="Y"&universe_df$lethal_mouse!=universe_df$lethal_het_mouse),]
naughties <- naughties[,c(2,33,34,43,51)]
colnames(naughties)[3:5] <- c("lethal_hom_mgi","lethal_hom_impc","lethal_het_mgi")

a<-lapply(naughties$gene,function(x) length(which(mgi_df$human_symbol==x)))
b<-lapply(naughties$gene,function(x) length(which(het$human_symbol==x)))
nos<-lapply(1:8, function(x) ifelse(a[[x]]>=b[[x]],a[[x]],b[[x]]))

allele_phenotype_check <- data.frame(gene = rep(naughties$gene,nos), MGI_ID = rep(naughties$MGI_ID,nos),
                                     lethal_hom_mgi = rep(naughties$lethal_hom_mgi,nos),lethal_het_mgi = rep(naughties$lethal_het_mgi,nos),
                                     hom_allele_compositions = rep(NA,sum(unlist(nos))),
                                     hom_allele_info = rep(NA,sum(unlist(nos))),
                                     hom_allele_type = rep(NA,sum(unlist(nos))),
                                     hom_allele_ID = rep(NA,sum(unlist(nos))),
                                     MP_ID_hom = rep(NA,sum(unlist(nos))),
                                     MP_DEF_hom = rep(NA,sum(unlist(nos))),
                                     het_allele_compositions = rep(NA,sum(unlist(nos))),
                                     het_allele_info = rep(NA,sum(unlist(nos))),
                                     het_allele_type = rep(NA,sum(unlist(nos))),
                                     het_allele_ID = rep(NA,sum(unlist(nos))),
                                     MP_ID_het = rep(NA,sum(unlist(nos))),
                                     MP_DEF_het = rep(NA,sum(unlist(nos)))
                                     )


homallcom<-lapply(naughties$gene, function(x) mgi_df$Allelic.Composition[which(mgi_df$human_symbol==x)])
allele_phenotype_check$hom_allele_compositions<-unlist(lapply(1:length(homallcom),function(x) {length(homallcom[[x]]) <- nos[[x]];homallcom[[x]]}))

homallinfo<-lapply(naughties$gene, function(x) mgi_df$allele_info[which(mgi_df$human_symbol==x)])
allele_phenotype_check$hom_allele_info<-unlist(lapply(1:length(homallinfo),function(x) {length(homallinfo[[x]]) <- nos[[x]];homallinfo[[x]]}))

homalltype<-lapply(naughties$gene, function(x) mgi_df$allele_type[which(mgi_df$human_symbol==x)])
allele_phenotype_check$hom_allele_type<-unlist(lapply(1:length(homalltype),function(x) {length(homalltype[[x]]) <- nos[[x]];homalltype[[x]]}))

homallID<-lapply(naughties$gene, function(x) mgi_df$allele_ID[which(mgi_df$human_symbol==x)])
allele_phenotype_check$hom_allele_ID<-unlist(lapply(1:length(homallID),function(x) {length(homallID[[x]]) <- nos[[x]];homallID[[x]]}))

homallmp<-lapply(naughties$gene, function(x) mgi_df$MP_ID[which(mgi_df$human_symbol==x)])
allele_phenotype_check$MP_ID_hom<-unlist(lapply(1:length(homallmp),function(x) {length(homallmp[[x]]) <- nos[[x]];homallmp[[x]]}))

homallmpdef<-lapply(naughties$gene, function(x) mgi_df$MP_DEF[which(mgi_df$human_symbol==x)])
allele_phenotype_check$MP_DEF_hom<-unlist(lapply(1:length(homallmpdef),function(x) {length(homallmpdef[[x]]) <- nos[[x]];homallmpdef[[x]]}))


hetallcom<-lapply(naughties$gene, function(x) het$Allelic.Composition[which(het$human_symbol==x)])
allele_phenotype_check$het_allele_compositions<-unlist(lapply(1:length(hetallcom),function(x) {length(hetallcom[[x]]) <- nos[[x]];hetallcom[[x]]}))

hetallinfo<-lapply(naughties$gene, function(x) het$allele_info[which(het$human_symbol==x)])
allele_phenotype_check$het_allele_info<-unlist(lapply(1:length(hetallinfo),function(x) {length(hetallinfo[[x]]) <- nos[[x]];hetallinfo[[x]]}))

hetalltype<-lapply(naughties$gene, function(x) het$allele_type[which(het$human_symbol==x)])
allele_phenotype_check$het_allele_type<-unlist(lapply(1:length(hetalltype),function(x) {length(hetalltype[[x]]) <- nos[[x]];hetalltype[[x]]}))

hetallID<-lapply(naughties$gene, function(x) het$allele_ID[which(het$human_symbol==x)])
allele_phenotype_check$het_allele_ID<-unlist(lapply(1:length(hetallID),function(x) {length(hetallID[[x]]) <- nos[[x]];hetallID[[x]]}))

hetallmp<-lapply(naughties$gene, function(x) het$MP_ID[which(het$human_symbol==x)])
allele_phenotype_check$MP_ID_het<-unlist(lapply(1:length(hetallmp),function(x) {length(hetallmp[[x]]) <- nos[[x]];hetallmp[[x]]}))

hetallmpdef<-lapply(naughties$gene, function(x) het$MP_DEF[which(het$human_symbol==x)])
allele_phenotype_check$MP_DEF_het<-unlist(lapply(1:length(hetallmpdef),function(x) {length(hetallmpdef[[x]]) <- nos[[x]];hetallmpdef[[x]]}))


write.xlsx(allele_phenotype_check, "allele_phenotype_check.xlsx")
