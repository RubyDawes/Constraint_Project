#a script to parse MGI_PhenoGenoMP to extract all genes with lethal phenotypes, write to xlsx
#load in all genes and phenotypes

mgi_df <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenoGenoMP.xlsx")
mgi_df <- mgi_df[,c(1,2,6,4)]
mgi_df<- mgi_df[which(!grepl(pattern = "\\,", x = mgi_df$MGI_ID, ignore.case = FALSE)),]

allele <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenotypicAllele.xlsx")
mgi_df$allele_info <- vlookup(mgi_df$`Allele.Symbol(s)`,allele,result_column="Allele.Attribute",lookup_column="Allele.Symbol")
mgi_df <- mgi_df[which(grepl(pattern="Null/knockout",x=mgi_df$allele_info,ignore.case=FALSE)|grepl(pattern="Hypomorph",x=mgi_df$allele_info,ignore.case=TRUE)),]

MGI_ID <- unique(mgi_df$MGI_ID)
mgi <- data.frame(MGI_ID)

#which genes have a lethal phenotype
mgi_lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
mgi_df$is_lethal <- ifelse(mgi_df$MP_ID%in%mgi_lethalphens$MP.id,"Y","N")

#getting all other phens
mgi_allphens <- read.xlsx("Gene_lists/MGI_MPs/VOC_MammalianPhenotype.xlsx",startRow=1,colNames=FALSE)
names(mgi_allphens)<- c("MP.id","phenotype","definition")

#high level phens
mgi_highphens <- read.xlsx("Gene_lists/MGI_MPs/HMD_HumanPhenotype.xlsx")
mgi_highphens <- mgi_highphens[,c(6,7)]
names(mgi_highphens) <- c("MGI.id","phenotype")
#mp ids
mgi$all_MP_ID <- lapply(mgi$MGI_ID, function(x) unique(mgi_df$MP_ID[which(mgi_df$MGI_ID==x)]))
mgi$lethal_MP_ID <- lapply(mgi$MGI_ID, function(x) unique(mgi_df$MP_ID[which(mgi_df$MGI_ID==x&mgi_df$is_lethal=="Y")]))
mgi$all_MP_phen <- lapply(mgi$all_MP_ID, function(x) mgi_allphens$phenotype[which(mgi_allphens$MP.id%in%x)])
mgi$lethal_MP_phen <- lapply(mgi$lethal_MP_ID, function(x) mgi_lethalphens$phenotype[which(mgi_lethalphens$MP.id%in%x)])

mgi$high_MP_ID <- strsplit(vlookup(mgi$MGI_ID,mgi_highphens,result_column="phenotype",lookup_column="MGI.id"),split=" ")
mgi$high_MP_phen <- lapply(mgi$high_MP_ID,function(x) unique(mgi_allphens$phenotype[which(mgi_allphens$MP.id%in%x)]))

#allele info
mgi$allele_info <- lapply(mgi$MGI_ID, function(x) unique(mgi_df$allele_info[which(mgi_df$MGI_ID==x&mgi_df$is_lethal=="Y")]))

#is_lethal
mgi$is_lethal <- ifelse(lapply(mgi$lethal_MP_ID,length)==0,"N","Y")

# add on gene names and human ortholog names
mgi_names <- read.xlsx("Gene_lists/MGI_MPs/HMD_HumanPhenotype.xlsx")
mgi_names$MGI.Marker.Accession.ID <- gsub(" ", "", mgi_names$MGI.Marker.Accession.ID, fixed = TRUE)
mgi$mouse_symbol <- vlookup(mgi$MGI_ID,mgi_names,result_column="Mouse.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi$human_symbol <- vlookup(mgi$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi$human_symbol <- checkGeneSymbols(mgi$human_symbol,unmapped.as.na=FALSE,hgnc.table=hgnc.table)[[3]]

mgi <- mgi[-which(is.na(mgi$mouse_symbol)),]  

rm(allele,mgi_names,mgi_df,mgi_lethalphens,MGI_ID)
save(mgi, file="output/Data/mgi.rda", compress="bzip2")
write.xlsx(mgi,"output/spreadsheets/MGI_genes_with_phenotypes.xlsx")
