# venn diagram showing overlap of wang, blomen + hart essentialomes
# venn diagram wang blomen hart overall
overrideTriple<- "ye"

# which are reproducibly in essentialome? (of genes common to 3 papers)
common_genes <- intersect(intersect(wang$Gene,blomen$Gene),hart$Gene)
wang_essentialomec <- intersect(common_genes,wang_essentialome)
blomen_essentialomec <- intersect(common_genes,blomen_essentialome)
hart_essentialomec <- intersect(common_genes,hart_essentialome)
#jpeg('overlaps_essentialome_common_BWH.jpg')
draw.triple.venn(area1 = length(wang_essentialomec), area2 = length(blomen_essentialomec), area3 = length(hart_essentialomec), n12 = length(intersect(wang_essentialomec,blomen_essentialomec)), n23 = length(intersect(blomen_essentialomec,hart_essentialomec)), n13 = length(intersect(wang_essentialomec,hart_essentialomec)), 
                 n123 = length(intersect(intersect(wang_essentialomec,blomen_essentialomec),hart_essentialomec)), category = c("Wang", "Blomen", "Hart"), lty = "blank", 
                 fill = c("dodgerblue3","brown3","slategray"), euler.d = TRUE, scaled = TRUE)
dev.off()
