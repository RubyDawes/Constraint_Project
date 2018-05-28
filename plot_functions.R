#theme for pie charts

blank_theme <- function() {
  theme_minimal(base_size=12,base_family="Avenir") %+replace%
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold"),
      legend.title=element_blank()
    )
}

density_theme <- function() {
  theme_bw(base_size=12,base_family="Avenir") %+replace%
    theme(
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid=element_blank(),
      plot.title=element_text(size=14, face="bold"),
      legend.title=element_blank()
    ) 
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none"))),
    legend,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}
