#Table One####
#sheet one- gene discovery informatics toolkit

OUT <- createWorkbook()
addWorksheet(OUT,"LEGEND")
addWorksheet(OUT,"S1_GD_informatics_toolkit")
addWorksheet(OUT,"S2_cell_ess_mouse_nonlethal")
addWorksheet(OUT,"S3_candidate_lethal_genes")

legend1<-read.table("output/spreadsheets/Table_One_Legend.txt",sep = "\t",fill = TRUE,quote = "")
colnames(legend1) <- legend1[1,1]
legend1 <- legend1[-1,]
legendstyle <- createStyle(fgFill = "#FFFF00", border = "TopBottom", textDecoration = "bold")
legendstyle1 <- createStyle(textDecoration = "bold")
addStyle(OUT,sheet = "LEGEND",style = legendstyle, rows = c(1,64,71),cols = c(1,2,3),gridExpand = TRUE)
addStyle(OUT,sheet = "LEGEND",style = legendstyle1, rows = c(2,65),cols = c(1,2,3),gridExpand = TRUE)

writeData(OUT,sheet = "LEGEND", x = legend1)
          

#sheet one- GD informatics toolkit
load("output/Data/universe_df.rda")
universe_df <- universe_df[-c(13:17)]
colnames(universe_df)[13]<-"gnomAD"

universestyle_gene <- createStyle(fgFill = "#CDCDCD", border = "TopBottom", textDecoration = "bold")
addStyle(OUT,sheet = "S1_GD_informatics_toolkit",style = universestyle_gene, rows = c(1),cols = c(1:4),gridExpand = TRUE)

universestyle_omim <- createStyle(fgFill = "#FFDA6B", border = "TopBottom", textDecoration = "bold")
addStyle(OUT,sheet = "S1_GD_informatics_toolkit",style = universestyle_omim, rows = c(1),cols = c(5:12),gridExpand = TRUE)

universestyle_constraint<- createStyle(fgFill = "#548AC5", border = "TopBottom", textDecoration = "bold")
addStyle(OUT,sheet = "S1_GD_informatics_toolkit",style = universestyle_constraint, rows = c(1),cols = c(13:26),gridExpand = TRUE)

universestyle_homko<- createStyle(fgFill = "#4EAC5B", border = "TopBottom", textDecoration = "bold")
addStyle(OUT,sheet = "S1_GD_informatics_toolkit",style = universestyle_homko, rows = c(1),cols = c(27:44),gridExpand = TRUE)

universestyle_hetko<- createStyle(fgFill = "#9FCD62", border = "TopBottom", textDecoration = "bold")
addStyle(OUT,sheet = "S1_GD_informatics_toolkit",style = universestyle_hetko, rows = c(1),cols = c(45:58),gridExpand = TRUE)

universestyle_cellko<- createStyle(fgFill = "#FFFF88", border = "TopBottom", textDecoration = "bold")
addStyle(OUT,sheet = "S1_GD_informatics_toolkit",style = universestyle_cellko, rows = c(1),cols = c(59:61),gridExpand = TRUE)

writeData(OUT,sheet = "S1_GD_informatics_toolkit", x = universe_df)

#sheet two- cell essential mouse non-lethal genes- info on phenotypes
cellnotmouse<-universe_df[which(universe_df$lethal_mouse=="N"&universe_df$cell_essential=="Y"),c("gene","MGI_ID","all_MP_phen","high_MP_phen")]
cellnotmouse$prem_death <- ifelse(grepl("premature death",cellnotmouse$all_MP_phen),"Y","N")
cellnotmouse$cellular_phenotype <- ifelse(grepl("cellular phenotype",cellnotmouse$high_MP_phen),"Y","N")
writeData(OUT,sheet = "S2_cell_ess_mouse_nonlethal", x = cellnotmouse)

#sheet three- 3,435 candidate lethal genes
candidates <- universe_df$gene[which(is.na(universe_df$omim)&(universe_df$lethal_mouse=="Y"|universe_df$lethal_het_mouse=="Y"|universe_df$cell_essential=="Y"))]
writeData(OUT,sheet = "S3_candidate_lethal_genes", x = candidates)

#saving table one
saveWorkbook(OUT,"output/spreadsheets/Table_One.xlsx")

#Table Two####
OUT <- createWorkbook()
addWorksheet(OUT,"LEGEND")
addWorksheet(OUT,"S1_lethal_search_terms")
addWorksheet(OUT,"S2_lethal_search_hit_details")
addWorksheet(OUT,"S3_cell_KO_tallies")
addWorksheet(OUT,"S4_lethal_MP_terms")

#legend
legend2<-read.table("output/spreadsheets/Table_Two_Legend.txt",sep = "\t",fill = TRUE,quote = "")
colnames(legend2) <- legend2[1,1]
legend2 <- legend2[-1,]
legendstyle <- createStyle(fgFill = "#FFFF00", border = "TopBottom", textDecoration = "bold")
legendstyle1 <- createStyle(textDecoration = "bold")
addStyle(OUT,sheet = "LEGEND",style = legendstyle, rows = c(1,2,15,30),cols = c(1,2,3),gridExpand = TRUE)
addStyle(OUT,sheet = "LEGEND",style = legendstyle1, rows = c(3,31),cols = c(1,2,3),gridExpand = TRUE)

writeData(OUT,sheet = "LEGEND", x = legend2)

#sheet one- search terms for human lethal genes OMIM API search
search_term<-'(("died as infants") OR ("born died shortly after"~2) OR ("death first year of life"~3) OR ("death infancy"~4) OR ("death infant"~4) OR ("death in utero"~2) OR ("death occurred hours"~4) OR ("death after birth"~3) OR ("delivered dead"~3) OR ("death within birth"~4) OR ("death months life"~4) OR ("died infants"~3) OR ("died days"~3) OR ("died before the age") OR ("died infancy"~3) OR ("died months life"~4) OR  ("dead months life"~4) OR ("death neonatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("death perinatal"~4) OR ("died perinatal"~4) OR ("dead perinatal"~4) OR ("lethal before birth"~4) OR ("lethal infantile"~4) OR ("lethal perinatal"~4) OR ("lethal prenatal"~4) OR ("lethal early life"~4) OR ("lethal first year of life"~3) OR ("lethal in utero"~2) OR ("die in utero"~2) OR ("lethal malformation"~3) OR ("lethal neonatal"~3) OR ("no development milestones"~2) OR ("severe lethal"~2) OR ("spontaneous abortion") OR ("the infant died") OR ("fatal infancy"~3) OR ("fetal demise"~3) OR ("embryonic lethal"~3) OR ("early lethal"~2) OR ("die after birth"~2) OR ("Died at weeks"~3) OR ("died at months"~4) OR ("die early childhood"~3) OR ("neonatally lethal"~2) OR ("lethal disorder"~4) OR ("died shortly after birth"~2) OR ("died hours after birth"~5) OR ("neonatal severe"~3) OR (lethal AND congenital) OR (lethal AND prenatal) OR (lethal AND perinatal) OR (lethal AND fetal) OR (lethal AND embryonic) OR (lethal AND neonatal) OR (lethal AND infantile) OR ("embryonically lethal"~3) OR ("prenatally lethal"~3) OR ("infants die"~3) OR ("die within the first year of life") OR ("fetal death") OR ("infants died"~3) OR ("died hours after delivery"~5) OR ("died 1 year"~3) OR ("died at months of age"~3) OR ("died in infancy"~4) OR ("stillborn fetus"~3) OR ("die first year of life"~3) OR ("died at months"~3) OR ("died days"~4) OR ("death at days"~3) OR ("death children"~5) OR ("death occurred month-old"~3) OR ("death early life"~3) OR ("death age months"~4) OR ("death by months"~5) OR ("died hours after birth"~7) OR ("died within months"~2) OR ("infant fatal"~3) OR ("infant died") OR ("early death") OR ("died early age"~3) OR ("died at hours"~2) OR ("died days of life"~3) OR ("death at months"~4))'
search_terms<-unlist(strsplit(search_term, "OR"))
search_terms <- unlist(lapply(search_terms,function(x) gsub('\\\"',"",x)))
writeData(OUT,sheet = "S1_lethal_search_terms", x = search_terms)

# sheet two- hits for human lethal genes OMIM API search

load("output/Data/human_lethal_hits.rda")
writeData(OUT,sheet = "S2_lethal_search_hit_details", x = results)


#sheet two- tallies for cell essential genes throughout 11 cell lines
load("output/Data/cell_KOs.rda")
write.xlsx(cell_KOs,"output/spreadsheets/Table_TwoS3.xlsx",sheetName = "S3_cell_KO_tallies",append=TRUE)
writeData(OUT,sheet = "S3_cell_KO_tallies", x = cell_KOs)

#sheet three- MP lethal terms
mgi_lethalphens <- read.xlsx("Gene_lists/MGI_MPs/MGI_lethal_phenotypes.xlsx")
writeData(OUT,sheet = "S4_lethal_MP_terms", x = mgi_lethalphens)

#saving table two
saveWorkbook(OUT,"output/spreadsheets/Table_Two.xlsx")
