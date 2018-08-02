universe_df$no_hpo_cats <- lengths(universe_df$hpo_slim)
ggplot(universe_df[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"),],aes(factor(cell_essential),no_hpo_cats))+
  geom_violin(aes(fill=factor(cell_essential)))

mean(universe_df$no_hpo_cats[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="Y")])
mean(universe_df$no_hpo_cats[which(universe_df$omim=="Y"&universe_df$lethal_mouse=="Y"&universe_df$cell_essential=="N")])


library(purrr)
library(readxl)

file <- "Gene_lists/HDN/HDN.xlsx"
sheets <- excel_sheets(file)
HDN <- map_df(sheets, ~ read_excel(file, sheet = .x,col_names=c("Disease ID","Disorder name","Gene symbols","OMIM ID","Chromosome","Class")))
rm(file,sheets)

HDN_gene_classes <- data.frame(class=unique(HDN$Class)[-c(9,23)])
HDN_gene_classes$genes <- lapply(HDN_gene_classes$class,function(x) 
  unique(unlist(strsplit(paste(
  unique(checkGeneSymbols(unlist(unique(strsplit(paste(HDN$`Gene symbols`[which(HDN$Class==x)],collapse=", "),", "))),unmapped.as.na=FALSE)[[3]])
  ,collapse=" /// ")," /// "))
  ))

save(HDN_gene_classes, file="output/Data/HDN_gene_classes.rda", compress="bzip2")
lengths(HDN_gene_classes$genes)
length(which(HDN_gene_classes$genes[[3]]%in%universe_df$gene))

universe_df$HDN_gene_class <-

which(universe_df$gene%in%HDN_gene_classes$genes[[1]])


