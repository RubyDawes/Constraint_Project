rec<-read.table("Gene_lists/BINGO/rec_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
rec$Genes.in.test.set <- lapply(rec$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))
rec$category <- 

dom<-read.table("Gene_lists/BINGO/dom_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
dom$Genes.in.test.set <- lapply(dom$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

candidates<-read.table("Gene_lists/BINGO/candidates_0.025.bgo",fill=TRUE,comment.char = "!",header=TRUE,sep="\t")
candidates$Genes.in.test.set <- lapply(candidates$Genes.in.test.set,function(x) unlist(strsplit(as.character(x),"\\|")))

rec$GO.ID <- lapply(rec$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
dom$GO.ID <- lapply(dom$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))
candidates$GO.ID <- lapply(candidates$GO.ID ,function(x) paste0("GO:",paste(rep("0",(7-nchar(x))),collapse=""),x))

goslim <- get_ontology("Gene_lists/GOSLIM/go.obo", propagate_relationships = "is_a",extract_tags = "minimal")

category <- unlist(lapply(goslim$ancestors,function(x) ifelse("GO:0003674"%in%x,"MF",
                                                                     ifelse("GO:0008150"%in%x,"BP",
                                                                            ifelse("GO:0005575"%in%x,"CC",NA)))))
catdf <- data.frame(id = names(category),cat = unname(category))
rm(category,goslim)

rec$category <- vlookup(unlist(rec$GO.ID),catdf,lookup_column = "id",result_column = "cat")
rec$category[which(is.na(rec$category))] <- c("BP","CC","CC")
dom$category <- vlookup(unlist(dom$GO.ID),catdf,lookup_column = "id",result_column = "cat")
dom$category[which(is.na(dom$category))] <- c("MF","CC","CC")
candidates$category <- vlookup(unlist(candidates$GO.ID),catdf,lookup_column = "id",result_column = "cat")
candidates$category[which(is.na(candidates$category))] <- c("CC","MF","BP")

rm(catdf)

rec <- rec[order(rec$category),]
dom <- dom[order(dom$category),]
candidates <- candidates[order(candidates$category),]


dom$redundancy<-rep(NA,length(dom$GO.ID))
dom$redundancy[which(dom$category=="BP")]<-lapply(which(dom$category=="BP"),function(y)
  dom$Description[which(unlist(lapply(which(dom$category=="BP"),function(x)
    length(which(dom$Genes.in.test.set[[y]]%in%dom$Genes.in.test.set[[x]]))/dom$x[[y]]*100
  ))>90)])

dom$redundancy[which(dom$category=="CC")]<-lapply(which(dom$category=="CC"),function(y)
  dom$Description[which(unlist(lapply(which(dom$category=="CC"),function(x)
    length(which(dom$Genes.in.test.set[[y]]%in%dom$Genes.in.test.set[[x]]))/dom$x[[y]]*100
  ))>90)+length(which(dom$category=="BP"))])

dom$redundancy[which(dom$category=="MF")]<-lapply(which(dom$category=="MF"),function(y)
  dom$Description[which(unlist(lapply(which(dom$category=="MF"),function(x)
    length(which(dom$Genes.in.test.set[[y]]%in%dom$Genes.in.test.set[[x]]))/dom$x[[y]]*100
  ))>90)+length(which(dom$category=="BP"))+length(which(dom$category=="CC"))])


rec$redundancy<-rep(NA,length(rec$GO.ID))
rec$redundancy[which(rec$category=="BP")]<-lapply(which(rec$category=="BP"),function(y)
  rec$Description[which(unlist(lapply(which(rec$category=="BP"),function(x)
    length(which(rec$Genes.in.test.set[[y]]%in%rec$Genes.in.test.set[[x]]))/rec$x[[y]]*100
  ))>90)])

rec$redundancy[which(rec$category=="CC")]<-lapply(which(rec$category=="CC"),function(y)
  rec$Description[which(unlist(lapply(which(rec$category=="CC"),function(x)
    length(which(rec$Genes.in.test.set[[y]]%in%rec$Genes.in.test.set[[x]]))/rec$x[[y]]*100
  ))>90)+length(which(rec$category=="BP"))])

rec$redundancy[which(rec$category=="MF")]<-lapply(which(rec$category=="MF"),function(y)
  rec$Description[which(unlist(lapply(which(rec$category=="MF"),function(x)
    length(which(rec$Genes.in.test.set[[y]]%in%rec$Genes.in.test.set[[x]]))/rec$x[[y]]*100
  ))>90)+length(which(rec$category=="BP"))+length(which(rec$category=="CC"))])

candidates$redundancy<-rep(NA,length(candidates$GO.ID))
candidates$redundancy[which(candidates$category=="BP")]<-lapply(which(candidates$category=="BP"),function(y)
  candidates$Description[which(unlist(lapply(which(candidates$category=="BP"),function(x)
    length(which(candidates$Genes.in.test.set[[y]]%in%candidates$Genes.in.test.set[[x]]))/candidates$x[[y]]*100
  ))>90)])

candidates$redundancy[which(candidates$category=="CC")]<-lapply(which(candidates$category=="CC"),function(y)
  candidates$Description[which(unlist(lapply(which(candidates$category=="CC"),function(x)
    length(which(candidates$Genes.in.test.set[[y]]%in%candidates$Genes.in.test.set[[x]]))/candidates$x[[y]]*100
  ))>90)+length(which(candidates$category=="BP"))])

candidates$redundancy[which(candidates$category=="MF")]<-lapply(which(candidates$category=="MF"),function(y)
  candidates$Description[which(unlist(lapply(which(candidates$category=="MF"),function(x)
    length(which(candidates$Genes.in.test.set[[y]]%in%candidates$Genes.in.test.set[[x]]))/candidates$x[[y]]*100
  ))>90)+length(which(candidates$category=="BP"))+length(which(candidates$category=="MF"))])



length(which(candidates$Genes.in.test.set[[which(candidates$Description=="primary metabolic process")]]%in%candidates$Genes.in.test.set[[which(candidates$Description=="metabolic process")]]))

#how many cytoplasm in all other CCs combined

length(which(rec$Genes.in.test.set[[which(rec$Description=="cytoplasm")]]%in%unlist(rec$Genes.in.test.set[which(!rec$Description=="cytoplasm"&rec$category=="CC")])))
length(rec$Genes.in.test.set[[which(rec$Description=="cytoplasm")]])


length(which(candidates$Genes.in.test.set[[which(candidates$Description=="cytoplasm")]]%in%unlist(candidates$Genes.in.test.set[which(!candidates$Description=="cytoplasm"&candidates$category=="CC")])))
length(candidates$Genes.in.test.set[[which(candidates$Description=="cytoplasm")]])

        
length(which(candidates$Genes.in.test.set[[which(candidates$Description=="nucleic acid binding")]]%in%unlist(candidates$Genes.in.test.set[which(!candidates$Description=="nucleic acid binding"&candidates$category=="MF")])))
length(candidates$Genes.in.test.set[[which(candidates$Description=="nucleic acid binding")]])
