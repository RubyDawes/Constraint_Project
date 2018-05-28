source("4.Cell_Knockouts/cell_essential_genes")
#seeing the proportion of genes in each set that are essential at different no-cell-lines-cutoffs
cell_line_cutoff_levels <- seq(0, 13, by=1)
cutoffs <- c()
for (i in 1:length(cell_line_cutoff_levels)) {cutoffs[i]= length(which(cell_KOs$no_hits>=cell_line_cutoff_levels[i]))}
cutoffs_perc <- c()
for (j in 1:length(cell_line_cutoff_levels)) {cutoffs_perc[j]=cutoffs[j]/length(cell_KOs$no_hits)*100}

plot(cell_line_cutoff_levels,cutoffs)
rm(i,j,cutoffs,cutoffs_perc)