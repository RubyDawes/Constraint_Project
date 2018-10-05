blank_theme <- function() {
  theme_minimal(base_size=12,base_family="Helvetica") %+replace%
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=12, face="bold"),
      legend.title=element_blank()
    )
}

scatter_theme <- function() {
  theme_bw(base_size=12,base_family="Helvetica") %+replace%
    theme(
      panel.grid=element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      plot.title=element_text(size=12, face="bold"),
      axis.title=element_text(size=12),
      axis.text=element_text(size=12)
    )
}

bar_theme <- function() {
  theme_bw(base_size=12,base_family="Helvetica") %+replace%
    theme(
      axis.title.x = element_blank(),
      panel.grid=element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(),
      legend.title=element_blank()
    )
}




