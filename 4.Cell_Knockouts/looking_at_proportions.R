
a<-lapply(seq(1, 11, by=1),function(x)length(which(universe_df$cell_essential_hits>=x&universe_df$omim=="Y"))/
            length(which(universe_df$cell_essential_hits>=x)))

b <- lapply(seq(1, 11, by=1),function(x)length(which(universe_df$cell_essential_hits>=x&universe_df$lethal_mouse=="Y"))/
              length(which(universe_df$cell_essential_hits>=x&!is.na(universe_df$lethal_mouse))))

