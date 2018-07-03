library("VennDiagram")
library("gridExtra")
load("output/Data/universe_df.rda")

#how many genes have cell + mouse phen data?
length(which(!is.na(universe_df$lethal_mouse)&!is.na(universe_df$cell_essential_hits)))
#how many genes are essential in either cell or mouse- 3094
length(which(universe_df[which(!is.na(universe_df$lethal_mouse)&!is.na(universe_df$cell_essential_hits)),]$lethal_mouse=="Y"|universe_df[which(!is.na(universe_df$lethal_mouse)&!is.na(universe_df$cell_essential_hits)),]$cell_essential_hits>=3))

#overlap between cell + mouse lethality when using cutoff of 3 cell lines
mouse=textGrob("Mouse Lethal", gp=gpar(fontsize=16,fontfamily="Helvetica"))
cell=textGrob("Cell Essential", gp=gpar(fontsize=16,fontfamily="Helvetica"))
title3=textGrob("Cutoff at 3 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a3=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=3))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<3))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=3)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_3celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a3), top=title3,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 2 cell lines
title2=textGrob("Cutoff at 2 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a2=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=2))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=2)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<2))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=2)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=2)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_2celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a2), top=title2,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 1 cell lines
title1=textGrob("Cutoff at 1 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a1=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=1))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=1)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<1))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=1)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=1)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_1celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a1), top=title1,left=mouse,right=cell)
dev.off()


#overlap between cell + mouse lethality when using cutoff of 4 cell lines
png("Sandra_Figures/Figs/cellmousevenn_4celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title4=textGrob("Cutoff at 4 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a4=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=4))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=4)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<4))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=4)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=4)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_4celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a4), top=title4,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 5 cell lines
png("Sandra_Figures/Figs/cellmousevenn_5celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title5=textGrob("Cutoff at 5 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a5=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=5))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=5)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<5))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=5)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=5)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_5celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a5), top=title5,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 6 cell lines
png("Sandra_Figures/Figs/cellmousevenn_6celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title6=textGrob("Cutoff at 6 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a6=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=6))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=6)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<6))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=6)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=6)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_6celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a6), top=title6,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 7 cell lines
png("Sandra_Figures/Figs/cellmousevenn_7celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title7=textGrob("Cutoff at 7 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a7=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=7))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=7)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<7))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=7)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=7)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_7celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a7), top=title7,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 8 cell lines
png("Sandra_Figures/Figs/cellmousevenn_8celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title8=textGrob("Cutoff at 8 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a8=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=8))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=8)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<8))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=8)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=8)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_8celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a8), top=title8,left=mouse,right=cell)
dev.off()


#overlap between cell + mouse lethality when using cutoff of 9 cell lines
png("Sandra_Figures/Figs/cellmousevenn_9celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title9=textGrob("Cutoff at 9 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a9=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=9))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=9)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<9))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=9)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=9)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_9celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a9), top=title9,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 10 cell lines
png("Sandra_Figures/Figs/cellmousevenn_10celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title10=textGrob("Cutoff at 10 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a10=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=10))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=10)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<10))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=10)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=10)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_10celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a10), top=title10,left=mouse,right=cell)
dev.off()


#overlap between cell + mouse lethality when using cutoff of 11 cell lines
png("Sandra_Figures/Figs/cellmousevenn_11celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title11=textGrob("Cutoff at 11 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a11=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=11))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=11)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<11))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=11)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=11)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_11celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a11), top=title11,left=mouse,right=cell)
dev.off()

#overlap between cell + mouse lethality when using cutoff of 12 cell lines
png("Sandra_Figures/Figs/cellmousevenn_12celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title12=textGrob("Cutoff at 12 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a12=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=12))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=12)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<12))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=12)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=12)), category = c("", ""), lty = "blank", 
                      fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_12celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a12), top=title12,left=mouse,right=cell)
dev.off()


#overlap between cell + mouse lethality when using cutoff of 13 cell lines
png("Sandra_Figures/Figs/cellmousevenn_13celllines.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
title13=textGrob("Cutoff at 13 cell lines", gp=gpar(fontsize=20,fontfamily="Helvetica"))
a13=draw.pairwise.venn(area1 = length(which(universe_df$lethal_mouse=="N"&universe_df$cell_essential_hits>=13))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=13)), area2 = length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits<13))+length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=13)), cross.area=length(which(universe_df$lethal_mouse=="Y"&universe_df$cell_essential_hits>=13)), category = c("", ""), lty = "blank", 
                       fill = c("indianred2", "steelblue2"), euler.d = TRUE, scaled = TRUE,cat.default.pos='outer',cex=c(1.5,1.5,1.5),fontfamily="Helvetica",alpha=0.8)
png("Sandra_Figures/Figs/cellmousevenn_13celllines.png",width=1300,height=1000,type="quartz",res=150,bg = "transparent")
grid.arrange(gTree(children=a13), top=title13,left=mouse,right=cell)
dev.off()


rm(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,cell,mouse,title1,title2,title3,title4,title5,title6,title7,title8,title9,title10,title11,title12,title13)
