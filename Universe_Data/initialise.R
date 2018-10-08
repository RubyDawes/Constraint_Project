hgnc<- read.csv("Gene_lists/Universe/gene_with_protein_product.txt",sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
hgnc$symbol<- checkGeneSymbols(hgnc$symbol)[[3]]
universe_df <- data.frame(hgnc$hgnc_id,hgnc$symbol,hgnc$name,hgnc$gene_family,hgnc$mgd_id)
names(universe_df) <- c("hgnc_id","gene","gene_name","gene_family","mgd_id")
rm(hgnc)
