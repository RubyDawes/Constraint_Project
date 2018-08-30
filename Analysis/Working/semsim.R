putdom <- unlist(universe_df$go_terms[
  which(is.na(universe_df$omim)&universe_df$lethal_mouse=="Y"&universe_df$constrained=="Y"&lengths(universe_df$go_terms)>0)])

dom<-unlist(universe_df$go_terms[
  which(universe_df$human_lethal=="Y"&(universe_df$Inheritance_pattern=="AD"|universe_df$Inheritance_pattern=="XLd")&lengths(universe_df$go_terms)>0)])


d <- godata('org.Hs.eg.db', ont="MF", computeIC=FALSE)
mgoSim(putdom,dom,semData=d, measure="Wang")