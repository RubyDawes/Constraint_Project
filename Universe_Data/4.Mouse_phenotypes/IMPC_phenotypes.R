#old data from 2016 release
#impc_old<-read.xlsx("Gene_lists/IMPC/IMPC_data.xlsx")
#impc_old$Gene_symbol<- checkGeneSymbols(impc_old$Gene_symbol,unmapped.as.na=FALSE)[[3]]

#putting in release 8.0

impc_all<- read.csv("Gene_lists/IMPC/ALL_genotype_phenotype.csv",header=TRUE, sep = ",",fill=TRUE,stringsAsFactors = FALSE)
impc_all<-impc_all[-which(impc_all$mp_term_id==""|impc_all$zygosity!="homozygote"),]

lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
impc_all$is_lethal <- ifelse(impc_all$mp_term_id%in%lethalphens$MP.id,"Y","N")

impc <- data.frame(mgi_id=unique(impc_all$marker_accession_id))

#mp ids
impc$all_MP_ID <- lapply(impc$mgi_id, function(x) unique(impc_all$mp_term_id[which(impc_all$marker_accession_id==x)]))
impc$all_MP_phen <- lapply(impc$mgi_id, function(x) unique(impc_all$mp_term_name[which(impc_all$marker_accession_id==x)]))
impc$lethal_MP_ID <- lapply(impc$mgi_id, function(x) unique(impc_all$mp_term_id[which(impc_all$marker_accession_id==x&impc_all$is_lethal=="Y")]))
impc$lethal_MP_phen <- lapply(impc$mgi_id, function(x) unique(impc_all$mp_term_name[which(impc_all$marker_accession_id==x&impc_all$is_lethal=="Y")]))

impc$is_lethal <- rep("N",length(impc$mgi_id))
impc$is_lethal[which(lengths(impc$lethal_MP_ID)>0)]<-"Y"

# add on gene names and human ortholog names
names <- read.csv("Gene_lists/Universe/gene_with_protein_product.txt",sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
impc$human_symbol <- vlookup(impc$mgi_id,names,result_column="symbol",lookup_column = "mgd_id")
impc$human_symbol <- checkGeneSymbols(impc$human_symbol,unmapped.as.na=FALSE)[[3]]
#removing pseudogenes/non-coding genes etc
impc <- impc[-which(is.na(impc$human_symbol)),]

save(impc, file="output/Data/impc.rda", compress="bzip2")

