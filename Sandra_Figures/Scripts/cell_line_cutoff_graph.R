#looking specifically at number of genes hat are essential in different numbers of cell lines throughout 3 papers (11 cell lines)
cell_line_cutoff_levels <- seq(1, 11, by=1)
cell_line_cutoff_numbers <- lapply(cell_line_cutoff_levels,function(x) length(which(universe_df$cell_essential_hits >= x)))
cutoffs <- data.frame(cell_line_cutoff_levels,unlist(cell_line_cutoff_numbers))
names(cutoffs)<-c("levels","numbers")

chosen <- subset(cutoffs, levels== "3")


cell_cutoff_plot <- ggplot(cutoffs,aes(x=levels,y=numbers))+
  geom_point(size=4,color="sandybrown")+scatter_theme()+
  labs(x="Number of cell lines",y="Number of essential genes")+
  scale_y_continuous(breaks = pretty(cutoffs$numbers, n = 10))+
  scale_x_continuous(breaks = pretty(cutoffs$levels, n = 10))


ggsave("Sandra_Figures/Figs/fig3b_cell_cutoffs.pdf",height=6.5, width=6.5, units='in')

