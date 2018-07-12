
proportion_cell_pli_overlap <- function(cell_essential_hits){
  return(length(which(universe_df$cell_essential_hits>=cell_essential_hits&universe_df$pLI>=0.9))/
           length(which(universe_df$cell_essential_hits>=cell_essential_hits&!is.na(universe_df$pLI))))
}
cell_pli_overlaps <- data.frame(number_cell_lines<-seq(1,11))
cell_pli_overlaps$perc_essential <-lapply(seq(1,11),proportion_cell_pli_overlap)

proportion_cell_mouse_overlap <- function(cell_essential_hits){
  return(length(which(universe_df$cell_essential_hits>=cell_essential_hits&universe_df$lethal_mouse=="Y"))/
           length(which(universe_df$cell_essential_hits>=cell_essential_hits&!is.na(universe_df$lethal_mouse))))
}

cell_mouse_overlaps <- data.frame(number_cell_lines<-seq(1,11))
cell_mouse_overlaps$perc_essential <-lapply(seq(1,11),proportion_cell_mouse_overlap)

proportion_cell_missense_overlap <- function(cell_essential_hits){
  return(length(which(universe_df$cell_essential_hits>=cell_essential_hits&universe_df$mis_z>=3.09))/
           length(which(universe_df$cell_essential_hits>=cell_essential_hits&!is.na(universe_df$mis_z))))
}
cell_mis_overlaps <- data.frame(number_cell_lines<-seq(1,11))
cell_mis_overlaps$perc_essential <-lapply(seq(1,11),proportion_cell_missense_overlap)

proportion_cell_omim_overlap <- function(cell_essential_hits){
  return(length(which(universe_df$cell_essential_hits>=cell_essential_hits&universe_df$omim=="Y"))/
           length(which(universe_df$cell_essential_hits>=cell_essential_hits)))
}
cell_omim_overlaps <- data.frame(number_cell_lines<-seq(1,11))
cell_omim_overlaps$perc_essential <-lapply(seq(1,11),proportion_cell_omim_overlap)


length(which(universe_df$cell_essential_hits>=3&universe_df$omim=="Y"))
length(which(universe_df$cell_essential_hits>=3&is.na(universe_df$omim)))
length(which(universe_df$cell_essential_hits>=3&universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AD"))


