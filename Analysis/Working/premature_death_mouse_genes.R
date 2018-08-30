
universe_df$prem_death<- rep(NA,length(universe_df$gene))
universe_df$prem_death[which(grepl("MP:0002083",universe_df$all_MP_ID))]<-"Y"

length(which(universe_df$lethal_mouse=="N"&universe_df$prem_death=="Y"))
#307 prem death genes which don't have any other lethal phenotype

length(which(universe_df$lethal_mouse=="N"&universe_df$prem_death=="Y"&universe_df$omim=="Y"))/length(which(universe_df$lethal_mouse=="N"&universe_df$prem_death=="Y"))
#125 are OMIM- 41%
length(which(universe_df$lethal_mouse=="N"&universe_df$prem_death=="Y"&is.na(universe_df$omim)))/length(which(universe_df$lethal_mouse=="N"&universe_df$prem_death=="Y"))
#182 aren't OMIM- 59%

length(which(universe_df$lethal_mouse=="Y"&universe_df$omim=="Y"))/length(which(universe_df$lethal_mouse=="Y"))
#36% of mouse-lethal genes are OMIM genes
length(which(universe_df$lethal_mouse=="N"&universe_df$omim=="Y"))/length(which(universe_df$lethal_mouse=="N"))
#20% of mouse-nonlethal genes are OMIM genes

#How many genes are just incomplete penetrance
incom <- read.xlsx("Gene_lists/MGI_MPs/incomplete_penetrance_phens.xlsx")
incom<-incom$MP.id

universe_df$incom<-lapply(universe_df$lethal_MP_ID,function(x) unique(x%in%incom))

length(which(universe_df$incom=="TRUE"))
length(which(universe_df$incom=="TRUE"))/length(which(universe_df$lethal_mouse=="Y"))
length(which(universe_df$incom=="TRUE"&universe_df$omim=="Y"))/length(which(universe_df$incom=="TRUE"))
length(which(universe_df$incom=="FALSE"&universe_df$omim=="Y"))/length(which(universe_df$incom=="FALSE"))
#784 genes are lethal with only incomplete penetrance
#24% of all mouse lethal genes
#- 31% of these are OMIM genes
#37% complete penetrance only lethal genes are OMIM genes
