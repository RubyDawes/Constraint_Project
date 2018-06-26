#CLEAR ALL VARIABLES

rm(list = ls())
library(HGNChelper)
library(openxlsx)
library(expss)
library(ggplot2)
library(grid)
library(png)
library(gridExtra)
library(scales)
library(VennDiagram)
library(kableExtra)
library(reshape2)
library(data.table)

load("output/Data/hgnc.table.rda")

hgnc<- read.csv("Gene_lists/Universe/gene_with_protein_product.txt",sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
hgnc$symbol<- checkGeneSymbols(hgnc$symbol,hgnc.table=hgnc.table)[[3]]
universe_df <- data.frame(hgnc$hgnc_id,hgnc$symbol,hgnc$name,hgnc$gene_family,hgnc$mgd_id)
names(universe_df) <- c("hgnc_id","gene","gene_name","gene_family","mgd_id")
universe <- hgnc$symbol
rm(hgnc)
