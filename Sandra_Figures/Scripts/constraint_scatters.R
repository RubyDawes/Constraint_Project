
nonomim <- ggplot(universe_df[which(is.na(universe_df$omim)&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
<<<<<<< HEAD
      geom_point(size=0.1)+scatter_theme()+
      geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1.5)+
=======
      geom_point(size=0.1,color="slategray")+scatter_theme()+
      geom_vline(xintercept=3.09,linetype="dashed",color="orange",size=1.5)+
>>>>>>> parent of a11dca3... aesthetic changes to scatter plots
      geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1.5)+
      labs(x="Missense Z Score",y="pLI score")+ggtitle("Constraint scores of non-OMIM genes\n")+
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01))+
      scale_x_continuous(expand = c(0, 0), limits = c(-9, 14))
png("Sandra_Figures/Figs/scatter_nonomim.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
nonomim
dev.off()



omim <- ggplot(universe_df[which(universe_df$omim=="Y"&!is.na(universe_df$exac)),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,color="slategray4")+scatter_theme()+
  geom_vline(xintercept=3.09,linetype="dashed",color="orange",size=1.5)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1.5)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle("Constraint scores of OMIM genes\n")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) +
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14)) 

png("Sandra_Figures/Figs/scatter_omim.png",width=1000,height=1000,type="quartz",res=150,bg = "transparent")
omim
dev.off()
