goa <- read.table("Gene_lists/Gene_Ontology/goa_human.gaf",header=FALSE,sep="\t",comment.char="!",fill=TRUE)
goa <- goa[-which(grepl("NOT",goa$V4)),]
goa$V3 <- checkGeneSymbols(goa$V3,unmapped.as.na = FALSE)[[3]]
universe_df$go_terms <- lapply(universe_df$gene, function(x) as.character(unique(goa$V5[which(goa$V3%in%x)])))
universe_df$go_terms[which(lengths(universe_df$go_terms)==0)] <- list(NULL)
rm(goa)
go <- get_ontology("Gene_lists/Gene_Ontology/go.obo", propagate_relationships = "is_a",extract_tags = "minimal")

universe_df$go_names <- lapply(universe_df$go_terms, function(x) lapply(x,function(y) ifelse(length(which(go$id==y))==0,"",go$name[[y]])))
rm(go)
go_annotations <- data.frame(universe_df$gene,I(universe_df$go_terms),I(universe_df$go_names))
save(go_annotations, file="output/Data/GO_annotations.rda", compress="bzip2")
rm(go_annotations)