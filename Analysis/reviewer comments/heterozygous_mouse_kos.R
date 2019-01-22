mgi_df <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenoGenoMP.xlsx")
mgi_df <- mgi_df[,c(1,2,6,4)]

#removing any rows with a comma in the MGI ID column - indicates a multi-gene knockout (we are only interested in single gene knockouts)
mgi_df<- mgi_df[which(!grepl(pattern = "\\,", x = mgi_df$MGI_ID, ignore.case = FALSE)),]

#homozygotes but with multiple alleles
homdifall <- mgi_df[which(grepl("\\|",mgi_df$`Allele.Symbol(s)`)==TRUE&grepl("\\+",mgi_df$`Allele.Symbol(s)`==FALSE)),]

# keeping only rows with data on heterozygous KO mice (by looking for a "|" in the allele syndrome column)

het <- mgi_df[which(grepl("\\+",mgi_df$`Allele.Symbol(s)`)==TRUE),]

#extracting only non-wildtype allele
het$allele <- strsplit(het$`Allele.Symbol(s)`,"\\|")
het$allele <- lapply(1:length(het$allele), function(x) het$allele[[x]][which(grepl("\\+",het$allele[[x]])==FALSE)])
het<- het[-which(lengths(het$allele)!=1),]
het$allele <- unlist(het$allele)

#reading in allele data so we can filter to only include null/knockout and hypomorph alleles. 
## In addition, this filters out heterozygotes
allele <- read.xlsx("Gene_lists/MGI_MPs/MGI_PhenotypicAllele.xlsx")
het$allele_info <- vlookup(het$allele,allele,result_column="Allele.Attribute",lookup_column="Allele.Symbol")
het <- het[which(grepl(pattern="Null/knockout",x=het$allele_info,ignore.case=FALSE)|grepl(pattern="Hypomorph",x=het$allele_info,ignore.case=TRUE)),]

#initialising df of data on heterozygous KO mouse phenotype from MGI
MGI_ID <- unique(het$MGI_ID)
mgi_het <- data.frame(MGI_ID)


#which genes have a lethal phenotype
mgi_lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
het$is_lethal <- ifelse(het$MP_ID%in%mgi_lethalphens$MP.id,"Y","N")

#getting all other phens
mgi_allphens <- read.xlsx("Gene_lists/MGI_MPs/VOC_MammalianPhenotype.xlsx",startRow=1,colNames=FALSE)
names(mgi_allphens)<- c("MP.id","phenotype","definition")

#mp ids
mgi_het$all_MP_ID <- lapply(mgi_het$MGI_ID, function(x) unique(het$MP_ID[which(het$MGI_ID==x)]))
mgi_het$lethal_MP_ID <- lapply(mgi_het$MGI_ID, function(x) unique(het$MP_ID[which(het$MGI_ID==x&het$is_lethal=="Y")]))

mgi_het$all_MP_phen <- lapply(mgi_het$all_MP_ID, function(x) mgi_allphens$phenotype[which(mgi_allphens$MP.id%in%x)])
mgi_het$lethal_MP_phen <- lapply(mgi_het$lethal_MP_ID, function(x) mgi_lethalphens$phenotype[which(mgi_lethalphens$MP.id%in%x)])


#allele info
mgi_het$allele_info <- lapply(mgi_het$MGI_ID, function(x) unique(het$allele_info[which(het$MGI_ID==x&het$is_lethal=="Y")]))

#is_lethal
mgi_het$is_lethal <- ifelse(lapply(mgi_het$lethal_MP_ID,length)==0,"N","Y")

# add on gene names and human ortholog names
mgi_names <- read.xlsx("Gene_lists/MGI_MPs/HMD_HumanPhenotype.xlsx")
mgi_names$MGI.Marker.Accession.ID <- gsub(" ", "", mgi_names$MGI.Marker.Accession.ID, fixed = TRUE)
mgi_het$mouse_symbol <- vlookup(mgi_het$MGI_ID,mgi_names,result_column="Mouse.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi_het$human_symbol <- vlookup(mgi_het$MGI_ID,mgi_names,result_column="Human.Marker.Symbol",lookup_column = "MGI.Marker.Accession.ID")
mgi_het$human_symbol <- checkGeneSymbols(mgi_het$human_symbol,unmapped.as.na=FALSE)[[3]]

mgi_het <- mgi_het[-which(is.na(mgi_het$mouse_symbol)),]  

rm(allele,het,mgi_allphens,mgi_df,mgi_highphens,mgi_lethalphens,mgi_names,MGI_ID)

#IMPC phenotypes
impc_all<- read.csv("Gene_lists/IMPC/ALL_genotype_phenotype.csv",header=TRUE, sep = ",",fill=TRUE,stringsAsFactors = FALSE)
impc_all<-impc_all[-which(impc_all$mp_term_id==""|impc_all$zygosity!="heterozygote"),]

#up to here
lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
impc_all$is_lethal <- ifelse(impc_all$mp_term_id%in%lethalphens$MP.id,"Y","N")

impc_het <- data.frame(mgi_id=unique(impc_all$marker_accession_id))

#mp ids
impc_het$all_MP_ID <- lapply(impc_het$mgi_id, function(x) unique(impc_all$mp_term_id[which(impc_all$marker_accession_id==x)]))
impc_het$all_MP_phen <- lapply(impc_het$mgi_id, function(x) unique(impc_all$mp_term_name[which(impc_all$marker_accession_id==x)]))
impc_het$lethal_MP_ID <- lapply(impc_het$mgi_id, function(x) unique(impc_all$mp_term_id[which(impc_all$marker_accession_id==x&impc_all$is_lethal=="Y")]))
impc_het$lethal_MP_phen <- lapply(impc_het$mgi_id, function(x) unique(impc_all$mp_term_name[which(impc_all$marker_accession_id==x&impc_all$is_lethal=="Y")]))

impc_het$is_lethal <- rep("N",length(impc_het$mgi_id))
impc_het$is_lethal[which(lengths(impc_het$lethal_MP_ID)>0)]<-"Y"

# add on gene names and human ortholog names
names <- read.csv("Gene_lists/Universe/gene_with_protein_product.txt",sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
impc_het$human_symbol <- vlookup(impc_het$mgi_id,names,result_column="symbol",lookup_column = "mgd_id")
impc_het$human_symbol <- checkGeneSymbols(impc_het$human_symbol,unmapped.as.na=FALSE)[[3]]
#removing pseudogenes/non-coding genes etc
impc_het <- impc_het[-which(is.na(impc_het$human_symbol)),]


rm(impc_all,names,lethalphens)

#appending heterozygous phenotypes from MGI and IMPC to universe_df
#MGI mouse knockouts- lethal or non-lethal in a heterozygous KO mouse
universe_df$lethal_het_MGI <- vlookup(universe_df$gene,mgi_het,result_column="is_lethal",lookup_column="human_symbol")
universe_df$lethal_het_MP_ID <- vlookup(universe_df$gene,mgi_het,result_column="lethal_MP_ID",lookup_column="human_symbol")
universe_df$lethal_het_MP_phen <- vlookup(universe_df$gene,mgi_het,result_column="lethal_MP_phen",lookup_column="human_symbol")
universe_df$all_het_MP_ID <- vlookup(universe_df$gene,mgi_het,result_column="all_MP_ID",lookup_column="human_symbol")
universe_df$all_het_MP_phen <- vlookup(universe_df$gene,mgi_het,result_column="all_MP_phen",lookup_column="human_symbol")
universe_df$het_allele_info <- vlookup(universe_df$gene,mgi_het,result_column="allele_info",lookup_column="human_symbol")

#IMPC phenotype data
universe_df$lethal_het_IMPC <- vlookup(universe_df$gene,impc_het,result_column="is_lethal",lookup_column="human_symbol")
universe_df$IMPC_het_all_MP_ID <- vlookup(universe_df$gene,impc_het,result_column="all_MP_ID",lookup_column="human_symbol")
universe_df$IMPC_het_all_MP_phen <- vlookup(universe_df$gene,impc_het,result_column="all_MP_phen",lookup_column="human_symbol")
universe_df$IMPC_het_lethal_MP_ID <- vlookup(universe_df$gene,impc_het,result_column="lethal_MP_ID",lookup_column="human_symbol")
universe_df$IMPC_het_lethal_MP_phen <- vlookup(universe_df$gene,impc_het,result_column="lethal_MP_phen",lookup_column="human_symbol")

universe_df$IMPC_het_ko <- rep(NA,length(universe_df$gene))
universe_df$IMPC_het_ko[which(lengths(universe_df$IMPC_all_MP_ID)>0)] <-  "Y"

universe_df$mouse_het_ko <- ifelse(!is.na(universe_df$lethal_het_MGI),"Y",ifelse(!is.na(universe_df$lethal_het_IMPC),"Y",NA))
universe_df$lethal_het_mouse <- rep(NA,length(universe_df$gene))
universe_df$lethal_het_mouse[which(universe_df$mouse_het_ko=="Y")] <- "N"
universe_df$lethal_het_mouse[which(universe_df$lethal_het_MGI=="Y"|universe_df$lethal_het_IMPC=="Y")] <- "Y"

rm(mgi_het,impc_het)

#are all the genes lethal in heterozygous KO mice, also lethal in homozygous KO mice?
lethal_homhet <- universe_df[which(universe_df$lethal_mouse=="Y"&universe_df$lethal_het_mouse=="Y"),]
lethal_hetonly <- universe_df[which(universe_df$lethal_mouse=="N"&universe_df$lethal_het_mouse=="Y"),]
lethal_hetnohom <- universe_df[which(is.na(universe_df$lethal_mouse)&universe_df$lethal_het_mouse=="Y"),]
lethal_homonly <- universe_df[which(universe_df$lethal_mouse=="Y"&universe_df$lethal_het_mouse=="N"),]
lethal_homnohet <- universe_df[which(universe_df$lethal_mouse=="Y"&is.na(universe_df$lethal_het_mouse)),]

#making venn diagram
###Figure 1C: Venn diagram showing overlap between mouse lethal genes ####
png('Analysis/reviewer comments/mouse_heterozygousKO_lethality_overlap.png',width=30,height=30,units="cm",res=1000)
draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="Y"&!is.na(universe_df$lethal_het_mouse))),
                   area2 = length(which(!is.na(universe_df$lethal_mouse)&universe_df$lethal_het_mouse=="Y")), 
                   cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$lethal_het_mouse=="Y")), 
                   category = c("homozygous KO mice", "heterozygous KO mice"), lty = "blank",fill = c("#ef7b0b", "#082e66"),alpha=0.6, euler.d = TRUE, 
                   scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica",cat.cex=(c(2,2)),
                   cat.fontfamily="Helvetica",cat.pos=c(220,140),cat.dist=c(.03,.03))

dev.off()



