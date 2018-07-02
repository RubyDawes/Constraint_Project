

length(which(universe_df$Inheritance_pattern=="AD"&universe_df$lethal_mouse=="Y"))
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$lethal_mouse=="N"))

length(which(universe_df$Inheritance_pattern=="AR"&universe_df$lethal_mouse=="Y"))
length(which(universe_df$Inheritance_pattern=="AR"&universe_df$lethal_mouse=="N"))

length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$lethal_mouse=="Y"))
length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$lethal_mouse=="N"))

length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="Y"))
length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="N"))

length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="Y"))
length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="N"))

