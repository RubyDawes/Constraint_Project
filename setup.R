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

load("hgnc.table.rda")

hgnc<- read.csv("Universe/Symbol_Check/gene_with_protein_product.txt",sep = "\t", comment.char = "#",stringsAsFactors = FALSE)
universe <- hgnc$symbol
rm(hgnc)
