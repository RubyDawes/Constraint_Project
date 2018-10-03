

length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))/
length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$pLI)))

length(which(universe_df$Inheritance_pattern=="AR"&universe_df$constrained=="N"))/
length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))

length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9))/
length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$pLI)))

length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="Y"))
length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="N"))

length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="Y"))
length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$lethal_mouse=="N"))

inh_lof <- data.frame(MT=c(length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$pLI<0.9))/length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))),
                      AR=c(length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac)))),
                      "ARAD"=c(length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac))),
                               length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac)))),
                      AD=c(length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac))),
                           length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9))/length(which(universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac)))),
                      XL=c(length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac))),
                           length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&universe_df$pLI<0.9))/length(which(grepl("XL",universe_df$Inheritance_pattern)&!grepl("MT",universe_df$Inheritance_pattern)&!is.na(universe_df$exac)))))



length(which(grepl(pattern="MT",universe_df$Inheritance_pattern)&universe_df$pLI>=0.9&universe_df$mis_z>=3.09))
length(which(universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09))
length(which(universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09))

length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z<3.09))
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z>=3.09))
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09))
length(which(universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9&universe_df$mis_z<3.09))


length(which(universe_df$omim=="Y"&universe_df$pLI>=0.9&universe_df$mis_z<3.09))
length(which(universe_df$omim=="Y"&universe_df$pLI<0.9&universe_df$mis_z>=3.09))
length(which(universe_df$omim=="Y"&universe_df$pLI>=0.9&universe_df$mis_z>=3.09))
length(which(universe_df$omim=="Y"&universe_df$pLI<0.9&universe_df$mis_z<3.09))
a3=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3))+
                        length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)),
                      area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3))+
                        length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), 
                      cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), 
                      category = c("", ""), lty = "blank",fill = c("sandybrown", "firebrick4"), euler.d = TRUE, 
                      scaled = TRUE,cat.default.pos='outer',cex=c(5,5,5),fontfamily="Helvetica")


length(which(is.na(universe_df$omim)&universe_df$pLI>=0.9&universe_df$mis_z<3.09))
length(which(is.na(universe_df$omim)&universe_df$pLI<0.9&universe_df$mis_z>=3.09))
length(which(is.na(universe_df$omim)&universe_df$pLI>=0.9&universe_df$mis_z>=3.09))
length(which(is.na(universe_df$omim)&universe_df$pLI<0.9&universe_df$mis_z<3.09))

length(which(universe_df$omim=="Y"&universe_df$mis_z>=3.09))/length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))
length(which(is.na(universe_df$omim)&universe_df$mis_z>=3.09))/length(which(is.na(universe_df$omim)&!is.na(universe_df$exac)))

length(which(universe_df$omim=="Y"&universe_df$pLI>=0.9))/length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))
length(which(is.na(universe_df$omim)&universe_df$pLI>=0.9))/length(which(is.na(universe_df$omim)&!is.na(universe_df$exac)))





#2d figure legend
length(which(universe_df$omim=="Y"&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/
  length(which(universe_df$omim=="Y"&!is.na(universe_df$exac)))

length(which(universe_df$omim=="Y"&universe_df$Inheritance_pattern=="AR"&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))/
  length(which(universe_df$omim=="Y"&universe_df$mis_z<3.09&universe_df$any_reg_constraint=="N"))



