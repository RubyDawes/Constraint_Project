source("Analysis/OMIM_API/API_search_functions.R")

search_term<-'("gain of function" OR "gain-of-function" OR "dominant negative" OR "dominant-negative")'
gof_results<-search_omim(search_term,"tx_molecular_genetics")
gof_results$omim_field_text <- lapply(gof_results$mim_numbers,function(x) retrieve_field(x,"molecularGenetics"))

gof_results$mim_entry_type<-vlookup(gof_results$mim_numbers,mapping,result_column = "entry_type",lookup_column="mim_number")
gof_results$prefix<-unlist(lapply(gof_results$mim_numbers, retrieve_prefix))
genemaps<-  lapply(gof_results$mim_numbers[which(gof_results$prefix=="#")],retrieve_genemap)

gene_mim_numbers <- c()
mapping_keys <- c()
inheritances <- c()
gene_symbols <- c()
for (i in 1:length(which(gof_results$prefix=="#"))) {
  if (typeof(genemaps[[i]])=="character"){
    gene_mim_numbers[i] <- NA
    mapping_keys[i]<-NA
    inheritances[i]<-NA
    gene_symbols[i]<-genemaps[[i]]
  }
  else {
    gene_mim_numbers[i]<-paste(as.character(genemaps[[i]][1]$gene_mim_numbers),collapse=", ")
    mapping_keys[i]<-paste(as.character(genemaps[[i]][2]$phenotype_mapping_key),collapse=", ")
    inheritances[i]<-paste(as.character(genemaps[[i]][3]$inheritance),collapse=", ")
    gene_symbols[i]<-paste(as.character(genemaps[[i]][4]$gene_symbol),collapse=", ")
  }}
rm(i)

gof_results$gene_mim_number[which(gof_results$prefix=="#")]<- gene_mim_numbers
gof_results$phenotype_mapping_keys[which(gof_results$prefix=="#")]<- mapping_keys
gof_results$inheritance[which(gof_results$prefix=="#")]<- inheritances
gof_results$gene_symbol[which(gof_results$prefix=="#")]<- gene_symbols
gof_results$gene_mim_number[which(gof_results$prefix=="*")]<-as.character(gof_results$mim_numbers[which(gof_results$prefix=="*")])
gof_results$gene_mim_number[which(gof_results$prefix=="%"|gof_results$mim_entry_type=="predominantly phenotypes")]<-"genetic basis unknown"
gof_results$gene_mim_number[which(is.na(gof_results$gene_mim_number))]<-"chromosomal abnormality"
gof_results$gene_mim_number<-lapply(gof_results$gene_mim_number,unique)
rm(genemaps,gene_mim_numbers,mapping_keys,inheritances,gene_symbols)

gof_results$gene_mim_number <- unlist(lapply(gof_results$gene_mim_number, function(x) unique(strsplit(x,", "))))

mim_ids<-unique(unlist(gof_results$gene_mim_number[which(gof_results$gene_mim_number!="chromosomal abnormality"&gof_results$gene_mim_number!="genetic basis unknown")]))
length(which(mim_ids%in%universe_df$mim_number))

gof_dn_genes <- universe_df[which(universe_df$mim_number%in%mim_ids),]


save(gof_dn_genes, file="output/Data/gof_dn_genes.rda", compress="bzip2")
save(gof_results, file="output/Data/gof_results",compress="bzip2")

#graphing constraint
scatter <- ggplot(gof_dn_genes[which(!is.na(gof_dn_genes$exac)&gof_dn_genes$human_lethal=="Y"),],aes(x=mis_z,y=pLI))+
  geom_point(size=0.1,alpha=0.3)+scatter_theme_goodsizes()+
  geom_vline(xintercept=3.09,linetype="dashed",color="lightsalmon2",size=1)+
  geom_hline(yintercept=0.9,linetype="dashed",color="indianred3",size=1)+
  labs(x="Missense Z Score",y="pLI score")+ggtitle(paste("GOF/DN \n n=",length(which(is.na(gof_dn_genes$omim)&!is.na(gof_dn_genes$exac))),"\n"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01))+
  scale_x_continuous(expand = c(0, 0), limits = c(-9, 14))
ggsave("Sandra_Figures/Figs/fig1ai.pdf",height=8.9, width=8.9, units='cm')

inh_mis <- data.frame(MT=c(length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$mis_z>=3.09))/length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac))),
                           length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$mis_z<3.09))/length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac)))),
                      AR=c(length(which(gof_dn_genes$Inheritance_pattern=="AR"&gof_dn_genes$mis_z>=3.09))/length(which(gof_dn_genes$Inheritance_pattern=="AR"&!is.na(gof_dn_genes$exac))),
                           length(which(gof_dn_genes$Inheritance_pattern=="AR"&gof_dn_genes$mis_z<3.09))/length(which(gof_dn_genes$Inheritance_pattern=="AR"&!is.na(gof_dn_genes$exac)))),
                      "ARAD"=c(length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&gof_dn_genes$mis_z>=3.09))/length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&!is.na(gof_dn_genes$exac))),
                               length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&gof_dn_genes$mis_z<3.09))/length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&!is.na(gof_dn_genes$exac)))),
                      AD=c(length(which(gof_dn_genes$Inheritance_pattern=="AD"&gof_dn_genes$mis_z>=3.09))/length(which(gof_dn_genes$Inheritance_pattern=="AD"&!is.na(gof_dn_genes$exac))),
                           length(which(gof_dn_genes$Inheritance_pattern=="AD"&gof_dn_genes$mis_z<3.09))/length(which(gof_dn_genes$Inheritance_pattern=="AD"&!is.na(gof_dn_genes$exac)))),
                      XL=c(length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$mis_z>=3.09))/length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac))),
                           length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$mis_z<3.09))/length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac)))))

inh_mism <- melt(inh_mis)
inh_mism$mis_constraint <- rep(c("Missense constraint     ", "No missense constraint     "), 5)

f <- ggplot(dat=inh_mism, aes(x=variable, y=value, fill=mis_constraint))
f<- f+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
f<- f+labs(y = "Proportion")+scale_fill_manual(values=c("lightsalmon2","steelblue3"))+theme(legend.position="bottom")
#ggsave("Sandra_Figures/Figs/fig1bi.pdf",height=7, width=7, units='cm')
rm(inh_mis,inh_mism,f)
inh_lof <- data.frame(MT=c(length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$pLI>=0.9))/length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac))),
                           length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$pLI<0.9))/length(which(grepl(pattern="MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac)))),
                      AR=c(length(which(gof_dn_genes$Inheritance_pattern=="AR"&gof_dn_genes$pLI>=0.9))/length(which(gof_dn_genes$Inheritance_pattern=="AR"&!is.na(gof_dn_genes$exac))),
                           length(which(gof_dn_genes$Inheritance_pattern=="AR"&gof_dn_genes$pLI<0.9))/length(which(gof_dn_genes$Inheritance_pattern=="AR"&!is.na(gof_dn_genes$exac)))),
                      "ARAD"=c(length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&gof_dn_genes$pLI>=0.9))/length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&!is.na(gof_dn_genes$exac))),
                               length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&gof_dn_genes$pLI<0.9))/length(which(gof_dn_genes$Inheritance_pattern=="AR,AD"&!is.na(gof_dn_genes$exac)))),
                      AD=c(length(which(gof_dn_genes$Inheritance_pattern=="AD"&gof_dn_genes$pLI>=0.9))/length(which(gof_dn_genes$Inheritance_pattern=="AD"&!is.na(gof_dn_genes$exac))),
                           length(which(gof_dn_genes$Inheritance_pattern=="AD"&gof_dn_genes$pLI<0.9))/length(which(gof_dn_genes$Inheritance_pattern=="AD"&!is.na(gof_dn_genes$exac)))),
                      XL=c(length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$pLI>=0.9))/length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac))),
                           length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&gof_dn_genes$pLI<0.9))/length(which(grepl("XL",gof_dn_genes$Inheritance_pattern)&!grepl("MT",gof_dn_genes$Inheritance_pattern)&!is.na(gof_dn_genes$exac)))))


inh_lofm <- melt(inh_lof)
inh_lofm$lof_constraint <- rep(c("LoF constraint     ", "no LoF constraint     "), 5)

e <- ggplot(dat=inh_lofm, aes(x=variable, y=value, fill=lof_constraint))
e<- e+geom_bar(width = 0.8, stat = "identity",color="slategray",position=position_fill(reverse = TRUE))+scale_y_continuous(expand = c(0, 0)) +bar_theme()
e<- e+labs(y = "Proportion")+scale_fill_manual(values=c("indianred3","steelblue3"))+theme(legend.position="bottom")


#ggsave("Sandra_Figures/Figs/fig1bii.pdf",height=7, width=7, units='cm')
rm(inh_lof,inh_lofm,e)




