# a script to read OMIM genemap2.txt file (downloaded from omim.org/downloads)
#, get gene names with well-established disease phenotype, get inheritance and 
# write to .xlsx

# Read OMIM data
omim <- read.csv(file = "Gene_lists/OMIM/genemap2.txt", sep = "\t", comment.char = "#",stringsAsFactors = FALSE)

#taking only needed columns
omim <- data.frame(omim$Mim.Number,omim$Approved.Symbol,omim$Gene.Name,omim$Phenotypes,omim$Mouse.Gene.Symbol.ID)
omim <- setNames(omim, c("Mim.Number","Approved.Symbol","Gene.Name","Phenotypes","Mouse.Gene.Symbol.ID"))

#changing types to character
omim$Phenotypes <- levels(omim$Phenotypes)[omim$Phenotypes]
omim$Approved.Symbol <- levels(omim$Approved.Symbol)[omim$Approved.Symbol]
omim$Mouse.Gene.Symbol.ID <- levels(omim$Mouse.Gene.Symbol.ID)[omim$Mouse.Gene.Symbol.ID]

# keep only rows containing phenotype mapping key (3) "molecular basis is known"
omim <- omim[which(grepl(pattern="\\(3)",x=omim$Phenotypes,ignore.case=FALSE)==TRUE),]

#look through phenotypes, filter out any nondiseases, susceptibilities, provisional links, somatic mutations
for (j in 1:length(omim$Phenotypes)) {
  omim$Phenotypes[j] <- strsplit(toString(omim$Phenotypes[j]),"; ")
}
for (i in 1:length(omim$Phenotypes)) {
  omim$Phenotypes[[i]]<- omim$Phenotypes[[i]][which(!grepl(pattern="\\[",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="\\{",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="\\?",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="somatic",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="Somatic",x=omim$Phenotypes[[i]],ignore.case=FALSE)&grepl(pattern="\\(3)",x=omim$Phenotypes[[i]],ignore.case=FALSE))]
}
omim <- omim[-which(lengths(omim$Phenotypes)==0),]

#mapping mim numbers to hgnc approved symbols
symbols <- read.csv(file = "Gene_lists/OMIM/mim2gene.txt", sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
symbols <- symbols[which(grepl(pattern="gene",symbols$MIM.Entry.Type..see.FAQ.1.3.at.https...omim.org.help.faq.,ignore.case=FALSE)),]
omim$Gene <- vlookup(omim$Mim.Number,symbols,result_column="Approved.Gene.Symbol..HGNC.",lookup_column="MIM.Number")
omim$Gene <- ifelse(omim$Gene=="",ifelse(omim$Mim.Number%in%symbols$MIM.Number,omim$Approved.Symbol,NA),omim$Gene)
omim<- omim[-which(is.na(omim$Gene)),]

#updating gene names
omim$Gene <- checkGeneSymbols(omim$Gene,unmapped.as.na=TRUE)[[3]]
omim <- omim[-which(is.na(omim$Gene)),]
omim$Gene[which(omim$Gene=="Sep-12")]<- "SEPT12"
omim$Gene[which(omim$Gene=="Sep-09")]<- "SEPT9"


#checking for duplicate genes, concatenating rows with duplicates, deleting duplicates
dups <- which(duplicated(omim$Gene)==TRUE)
for (p in 1:length(dups)) {
  genedups <- omim$Gene[dups[p]]
  index <- which(omim$Gene==genedups)
  omim$Mim.Number[index[1]]<- paste(omim$Mim.Number[index[1]],omim$Mim.Number[index[2]],sep=", ")
  omim$Phenotypes[index[1]]<- paste(omim$Phenotypes[index[1]],omim$Phenotypes[index[2]],sep="; ")
  if (length(omim$Mouse.Gene.Symbol.ID[index[which(nchar(omim$Mouse.Gene.Symbol.ID[index])>0)]])>0) {
    omim$Mouse.Gene.Symbol.ID[index[1]]<- omim$Mouse.Gene.Symbol.ID[index[which(nchar(omim$Mouse.Gene.Symbol.ID[index])>0)]]
  }
}
if (length(dups)>0){
  omim <- omim[-dups,]
}

#removing non protein-coding genes (noncoding RNAs, complex loci, immunoglobulins etc)
omim <- omim[-which(!omim$Gene%in%universe&!grepl(pattern="///",omim$Gene)),]

#getting inheritances
# Extend OMIM data
omim$MT <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = ", Mitochondrial", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$AR <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "Autosomal recessive", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$AD <- ifelse(omim$Phenotypes == "", NA,ifelse(grepl(pattern = "Autosomal dominant", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$XLd <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "X-linked dominant", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$XLr <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "X-linked recessive", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$IS <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "isolated cases", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))

omim$Inheritance_pattern <- rep(NA,length(omim$Gene))
omim$Inheritance_pattern[which(omim$MT=="Y")] <- "MT"
omim$Inheritance_pattern[which(omim$XLd == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$XLd == "Y")]),"XLd",paste(omim$Inheritance_pattern[which(omim$XLd == "Y")],"XLd",sep=","))
omim$Inheritance_pattern[which(omim$XLr == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$XLr == "Y")]),"XLr",paste(omim$Inheritance_pattern[which(omim$XLr == "Y")],"XLr",sep=","))
omim$Inheritance_pattern[which(omim$AR == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$AR == "Y")]),"AR",paste(omim$Inheritance_pattern[which(omim$AR == "Y")],"AR",sep=","))
omim$Inheritance_pattern[which(omim$AD == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$AD == "Y")]),"AD",paste(omim$Inheritance_pattern[which(omim$AD == "Y")],"AD",sep=","))
omim$Inheritance_pattern[which(omim$IS == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$IS == "Y")]),"IS",paste(omim$Inheritance_pattern[which(omim$IS == "Y")],"IS",sep=","))


rm(dups,genedups,i,index,j,p,symbols)
#taking only needed rows, saving to spreadsheet
omim <- omim[,c(1,4,5,6,13)]
write.xlsx(omim, file= "output/spreadsheets/omim_filtered_with_inheritances.xlsx", row.names = FALSE,col.names = TRUE)



