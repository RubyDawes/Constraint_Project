#theme for pie charts

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
bar_twosided_theme <- function() {
  theme_bw(base_size=20,base_family="Helvetica") %+replace%
    theme(
      panel.grid=element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      legend.title=element_blank(),
      legend.position="none"
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
scatter_theme_goodsizes <- function() {
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

bar_blank_theme <- function() {
  theme_bw(base_size=40,base_family="Helvetica") %+replace%
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      panel.grid=element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line = element_line(),
      plot.title=element_text(size=14, face="bold"),
      legend.title=element_blank()
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
bar_theme_go <- function() {
  theme_bw(base_size=20,base_family="Helvetica") %+replace%
    theme(
      panel.grid=element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(),
      legend.title=element_blank()
    )
}
bar_theme_or <- function() {
  theme_bw(base_size=12,base_family="Helvetica") %+replace%
    theme(
      panel.grid=element_blank(),
      panel.border = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.line.y = element_line(),
      legend.title=element_blank(),
      legend.position = "none"
    )
}





density_theme <- function() {
  theme_bw(base_size=12,base_family="Helvetica") %+replace%
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
