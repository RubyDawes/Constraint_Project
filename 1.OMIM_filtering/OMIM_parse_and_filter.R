# a script to read OMIM genemap2.txt file (downloaded from omim.org/downloads)
#, get gene names with well-established disease phenotype, get inheritance and 
# write to .xlsx

# Read OMIM data
omim <- read.csv(file = "Gene_lists/OMIM/genemap2.txt", sep = "\t", comment.char = "#",stringsAsFactors = FALSE)

#taking only needed columns
omim <- data.frame(omim$Mim.Number,omim$Approved.Symbol,omim$Phenotypes,omim$Mouse.Gene.Symbol.ID)
omim <- setNames(omim, c("Mim.Number","Approved.Symbol","Phenotypes","Mouse.Gene.Symbol.ID"))

#changing types to character
omim$Phenotypes <- levels(omim$Phenotypes)[omim$Phenotypes]
omim$Mouse.Gene.Symbol.ID <- levels(omim$Mouse.Gene.Symbol.ID)[omim$Mouse.Gene.Symbol.ID]

# keep only rows containing phenotype mapping key (3) "molecular basis is known"
omim <- omim[which(grepl(pattern="\\(3)",x=omim$Phenotypes,ignore.case=FALSE)==TRUE),]

#updating gene names
omim$Approved.Symbol <- checkGeneSymbols(omim$Approved.Symbol,unmapped.as.na=FALSE)[[3]]

#checking for duplicate genes, concatenating rows with duplicates, deleting duplicates
dups <- which(duplicated(omim$Approved.Symbol)==TRUE)
for (p in 1:length(dups)) {
  genedups <- omim$Approved.Symbol[dups[p]]
  index <- which(omim$Approved.Symbol==genedups)
  omim$Mim.Number[index[1]]<- paste(omim$Mim.Number[index[1]],omim$Mim.Number[index[2]],sep=", ")
  omim$Phenotypes[index[1]]<- paste(omim$Phenotypes[index[1]],omim$Phenotypes[index[2]],sep="; ")
  if (length(omim$Mouse.Gene.Symbol.ID[index[which(nchar(omim$Mouse.Gene.Symbol.ID[index])>0)]])>0) {
    omim$Mouse.Gene.Symbol.ID[index[1]]<- omim$Mouse.Gene.Symbol.ID[index[which(nchar(omim$Mouse.Gene.Symbol.ID[index])>0)]]
  }
}
omim <- omim[-dups,]

#look through phenotypes, filter out any nondiseases, susceptibilities, provisional links, somatic mutations
for (j in 1:length(omim$Phenotypes)) {
  omim$Phenotypes[j] <- strsplit(toString(omim$Phenotypes[j]),"; ")
}
for (i in 1:length(omim$Phenotypes)) {
  omim$Phenotypes[[i]]<- omim$Phenotypes[[i]][which(!grepl(pattern="\\[",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="\\{",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="\\?",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="somatic",x=omim$Phenotypes[[i]],ignore.case=FALSE)&!grepl(pattern="Somatic",x=omim$Phenotypes[[i]],ignore.case=FALSE)&grepl(pattern="\\(3)",x=omim$Phenotypes[[i]],ignore.case=FALSE))]
}
omim <- omim[-which(lengths(omim$Phenotypes)==0),]

#getting inheritances
# Extend OMIM data
omim$MT <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = ", Mitochondrial", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$AR <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "Autosomal recessive", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$AD <- ifelse(omim$Phenotypes == "", NA,ifelse(grepl(pattern = "Autosomal dominant", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$XLd <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "X-linked dominant", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$XLr <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "X-linked recessive", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))
omim$IS <- ifelse(omim$Phenotypes == "",NA,ifelse(grepl(pattern = "isolated cases", x = omim$Phenotypes, ignore.case = FALSE),"Y","N"))

omim$Inheritance_pattern <- rep(NA,length(omim$Approved.Symbol))
omim$Inheritance_pattern[which(omim$MT=="Y")] <- "MT"
omim$Inheritance_pattern[which(omim$XLd == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$XLd == "Y")]),"XLd",paste(omim$Inheritance_pattern[which(omim$XLd == "Y")],"XLd",sep=","))
omim$Inheritance_pattern[which(omim$XLr == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$XLr == "Y")]),"XLr",paste(omim$Inheritance_pattern[which(omim$XLr == "Y")],"XLr",sep=","))
omim$Inheritance_pattern[which(omim$AR == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$AR == "Y")]),"AR",paste(omim$Inheritance_pattern[which(omim$AR == "Y")],"AR",sep=","))
omim$Inheritance_pattern[which(omim$AD == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$AD == "Y")]),"AD",paste(omim$Inheritance_pattern[which(omim$AD == "Y")],"AD",sep=","))
omim$Inheritance_pattern[which(omim$IS == "Y")]<- ifelse(is.na(omim$Inheritance_pattern[which(omim$IS == "Y")]),"IS",paste(omim$Inheritance_pattern[which(omim$IS == "Y")],"IS",sep=","))


names(omim)[names(omim) == 'Approved.Symbol'] <- 'gene'

rm(dups,genedups,i,index,j,p)
#taking only needed rows, saving to spreadsheet
omim <- omim[,c(1,2,3,4,11)]
write.xlsx(omim, file= "output/spreadsheets/omim_filtered_with_inheritances.xlsx", row.names = FALSE,col.names = TRUE)



