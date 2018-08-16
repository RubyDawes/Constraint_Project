omim_constraint_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$mis_z>=3.09)),
                             length(which(universe_df$omim=="Y"&universe_df$mis_z<3.09)),
                             length(which(is.na(universe_df$omim)&universe_df$mis_z>=3.09)),
                             length(which(is.na(universe_df$omim)&universe_df$mis_z<3.09)))

omim_constraint_OR<- OR_test(length(which(universe_df$omim=="Y"&universe_df$pLI>=0.9)),
                             length(which(universe_df$omim=="Y"&universe_df$pLI<0.9)),
                             length(which(is.na(universe_df$omim)&universe_df$pLI>=0.9)),
                             length(which(is.na(universe_df$omim)&universe_df$pLI<0.9)))
