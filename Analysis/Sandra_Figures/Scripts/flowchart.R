
#disease involvement row
length(which(universe_df$omim=="Y"))
length(which(is.na(universe_df$omim)))

#no constraint data
length(which(universe_df$omim=="Y"&is.na(universe_df$exac)))
length(which(is.na(universe_df$omim)&is.na(universe_df$exac)))

#ExAC constraint row
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"))
length(which(universe_df$omim=="Y"&universe_df$constrained=="N"))
length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"))
length(which(is.na(universe_df$omim)&universe_df$constrained=="N"))

#no mouse data
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&is.na(universe_df$mouse_ko)))
length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&is.na(universe_df$mouse_ko)))
length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&is.na(universe_df$mouse_ko)))
length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&is.na(universe_df$mouse_ko)))

#has mouse data
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&!is.na(universe_df$mouse_ko)))
length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&!is.na(universe_df$mouse_ko)))
length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&!is.na(universe_df$mouse_ko)))
length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&!is.na(universe_df$mouse_ko)))

#mouse KO phenotype
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y"))
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="N"))

length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&universe_df$lethal_mouse=="Y"))
length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&universe_df$lethal_mouse=="N"))

length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y"))
length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&universe_df$lethal_mouse=="N"))

length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&universe_df$lethal_mouse=="Y"))
length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&universe_df$lethal_mouse=="N"))

#percentages
length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$constrained=="Y"&!is.na(universe_df$mouse_ko)))

length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&universe_df$lethal_mouse=="Y"))/length(which(universe_df$omim=="Y"&universe_df$constrained=="N"&!is.na(universe_df$mouse_ko)))

length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$constrained=="Y"&!is.na(universe_df$mouse_ko)))

length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&universe_df$lethal_mouse=="Y"))/length(which(is.na(universe_df$omim)&universe_df$constrained=="N"&!is.na(universe_df$mouse_ko)))
