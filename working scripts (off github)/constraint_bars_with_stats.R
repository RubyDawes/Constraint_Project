source("2.ExAC_constraint/exac_constraint.R")

omim_exac <- exac[which(exac$disease=="Y"),]
nonomim_exac <- exac[which(exac$disease=="N"),]

mis <- data.frame(OMIM=c(length(which(omim_exac$mis_z>=3.09))/length(omim_exac$gene),length(which(omim_exac$mis_z<3.09))/length(omim_exac$gene)),non_OMIM=c(length(which(nonomim_exac$mis_z>=3.09))/length(nonomim_exac$gene),length(which(nonomim_exac$mis_z<3.09))/length(nonomim_exac$gene)))



########## Data Loading
dat <- read.xlsx("Week 5 - Sandra Cooper/MGI_stats_data.xlsx",
                 startRow = 1, colNames = TRUE, rowNames = FALSE)
dat_table <- table(dat$phenotype.in.MGI, dat$OMIM.gene) colnames(dat_table) = c("Non-OMIM", "OMIM") rownames(dat_table) = c("Lethal", "Non-Lethal") kable(dat_table, format = "latex", booktabs=T)
########## Chi Squared Test chisq.test(dat_table, correct = FALSE)
dat_df <- data.frame(
  OMIM = c(dat_table[1,2]/sum(dat_table[,2]), dat_table[2,2]/sum(dat_table[,2])), non_OMIM = c(dat_table[1,1]/sum(dat_table[,1]), dat_table[2,1]/sum(dat_table[,1])))
mism <- melt(mis)
mism$mis_constraint <- rep(c("missense constraint", "no missense constraint"), 2)
########## Plotting Lethality Proportions vs OMIM
ggplot(dat = mism, aes(x = variable, y = value, fill = mis_constraint)) +
  geom_bar(stat = "identity", colour = "black") +
  geom_text(position = "stack", aes(label = round(value,2)), color = "black", vjust =2) + scale_fill_manual(values = alpha(c("#3B9AB2", "#FF9999"),0.75)) +
  theme_minimal() +
  theme(axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black")) + labs(y = "Proportion", x="", fill = "Missense Constraint")
