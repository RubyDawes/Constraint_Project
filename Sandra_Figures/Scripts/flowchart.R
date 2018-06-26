length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3)) #4644
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3)) #125
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3)) #2118
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)) #851

universe_df$constrained <- ifelse(universe_df$mis_z>=3.09|universe_df$pLI>=0.9,"Y","N")
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3&universe_df$constrained=="N"))#3768
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3&universe_df$constrained=="Y"))#800

length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3&universe_df$constrained=="N"))#91
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3&universe_df$constrained=="Y"))#32

length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3&universe_df$constrained=="N"))#1199
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3&universe_df$constrained=="Y")) #868

length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3&universe_df$constrained=="N")) #480
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3&universe_df$constrained=="Y")) #364





length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3&universe_df$constrained=="N"&universe_df$omim=="Y"))#766
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3&universe_df$constrained=="N"&is.na(universe_df$omim)))#3185

length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3&universe_df$constrained=="Y"&universe_df$omim=="Y"))#
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits<3&universe_df$constrained=="Y"&is.na(universe_df$omim)))#

length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3&universe_df$constrained=="N"&universe_df$omim=="Y"))#26
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3&universe_df$constrained=="N"&is.na(universe_df$omim)))#65

length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3&universe_df$constrained=="Y"&universe_df$omim=="Y"))#5
length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3&universe_df$constrained=="Y"&is.na(universe_df$omim)))#27

length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3&universe_df$constrained=="N"&universe_df$omim=="Y"))#502
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3&universe_df$constrained=="N"&is.na(universe_df$omim)))#697

length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3&universe_df$constrained=="Y"&universe_df$omim=="Y")) #330
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3&universe_df$constrained=="Y"&is.na(universe_df$omim))) #538

length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3&universe_df$constrained=="N"&universe_df$omim=="Y")) #148
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3&universe_df$constrained=="N"&is.na(universe_df$omim))) #332

universe_df$gene_name[which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3&universe_df$constrained=="Y"&universe_df$omim=="Y")] #98
length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3&universe_df$constrained=="Y"&is.na(universe_df$omim))) #260

library(DiagrammeR)


nodes <- create_nodes(nodes=c("a","b","c","d"),
                      label = FALSE,
                      type = "lower",
                      style = "filled",
                      color = "aqua",
                      shape = c("circle", "circle",
                                "rectangle", "rectangle"),
                      data = c(3.5, 2.6, 9.4, 2.7))
  
grViz("
digraph neato {
graph [layout = neato]

      node [shape = circle,
      style = filled,
      color = grey,
      label = '']
      
      node [fillcolor = red]
      genes
      
      node [fillcolor = green]
      b c
      
      node [fillcolor = orange]
      
      edge [color = grey]
      genes -> {b c}
      b -> {e f g h i j}
      c -> {k l m n o p}
  }   
      ")

universe_df$constrained <- ifelse(universe_df$mis_z>=3.09|universe_df$pLI>=0.9,"Y","N")
length(which(universe_df$omim=="Y"))
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"))
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y"))
length(which(universe_df$omim=="Y"&universe_df$constrained=="N"))
length(which(is.na(universe_df$omim)))
length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"))
length(which(is.na(universe_df$omim)&universe_df$constrained=="N"))




