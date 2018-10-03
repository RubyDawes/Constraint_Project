# omim genes

unique(universe_df$Inheritance_pattern)

mt_no<-length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
               universe_df$Inheritance_pattern=="MT,AR,AD")&!is.na(universe_df$exac)))

length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$mis_z>=3.09))/mt_no*100

length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$pLI>=0.9))/mt_no*100

ar_no<-length(which((universe_df$Inheritance_pattern=="AR"&!is.na(universe_df$exac))))

length(which((universe_df$Inheritance_pattern=="AR"&universe_df$mis_z>=3.09)))/ar_no*100

length(which((universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9)))/ar_no*100

ar_ad_no<-length(which((universe_df$Inheritance_pattern=="AR,AD"&!is.na(universe_df$exac))))

length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$mis_z>=3.09)))/ar_ad_no*100

length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9)))/ar_ad_no*100

ad_no<-length(which((universe_df$Inheritance_pattern=="AD"&!is.na(universe_df$exac))))

length(which((universe_df$Inheritance_pattern=="AD"&universe_df$mis_z>=3.09)))/ad_no*100

length(which((universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9)))/ad_no*100


xl_no<-length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&!is.na(universe_df$exac))))

length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$mis_z>=3.09)))/xl_no*100

length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$pLI>=0.9)))/xl_no*100

# omim nonlethal genes

unique(universe_df$Inheritance_pattern)

mt_no<-length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                       universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac)))

length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))/mt_no*100

length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))/mt_no*100

ar_no<-length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))

length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/ar_no*100

length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/ar_no*100

ar_ad_no<-length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))

length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/ar_ad_no*100

length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/ar_ad_no*100

ad_no<-length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))

length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/ad_no*100

length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/ad_no*100


xl_no<-length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&!is.na(universe_df$exac))))

length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09)))/xl_no*100

length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9)))/xl_no*100


 # lethal genes

unlist(unique(universe_df$lethal_inheritance))

mt_no<-length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&!is.na(universe_df$exac)))

length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z>=3.09))/mt_no*100

length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$pLI>=0.9))/mt_no*100

ar_no<-length(which((universe_df$lethal_inheritance=="AR"&!is.na(universe_df$exac))))

length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z>=3.09)))/ar_no*100

length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI>=0.9)))/ar_no*100

ar_ad_no<-length(which((universe_df$lethal_inheritance=="AR,AD"&!is.na(universe_df$exac))))

length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z>=3.09)))/ar_ad_no*100

length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI>=0.9)))/ar_ad_no*100

ad_no<-length(which((universe_df$lethal_inheritance=="AD"&!is.na(universe_df$exac))))

length(which((universe_df$lethal_inheritance=="AD"&universe_df$mis_z>=3.09)))/ad_no*100

length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9)))/ad_no*100


xl_no<-length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                        universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&!is.na(universe_df$exac))))

length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                 universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z>=3.09)))/xl_no*100

length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                 universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI>=0.9)))/xl_no*100


#statistical tests

OR_test <- function(a,b,c,d){
  matrix <- matrix(c(a,b,c,d),nrow=2,dimnames=list(metric1 = c("Y","N"),metric2 = c("Y","N")))
  test <- fisher.test(matrix,alternative = "two.sided")
  return(test)
}

# MT MISSENSE
OR_test(length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                                             universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$mis_z>=3.09))+1,
                             length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                                             universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$mis_z<3.09))-1,
                             length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z>=3.09))+1,
                             length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$mis_z<3.09))-1)
# MT LOF
OR_test(length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                        universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$pLI>=0.9))+1,
        length(which((universe_df$Inheritance_pattern=="MT,AR"|universe_df$Inheritance_pattern=="MT,XLd,AR"|
                        universe_df$Inheritance_pattern=="MT,AR,AD")&universe_df$pLI<0.9))-1,
        length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$pLI>=0.9))+1,
        length(which((universe_df$lethal_inheritance=="MT,AR"|universe_df$lethal_inheritance=="MT,XLd")&universe_df$pLI<0.9))-1)

# AR MISSENSE
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$mis_z<3.09))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z<3.09)))
)

# AR LOF
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI<0.9)))
)

# AR,AD MISSENSE
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$mis_z<3.09))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z<3.09)))
)

# AR,AD LOF
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI<0.9)))  
)

# AD MISSENSE
OR_test(
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI<0.9)))
)

# AD LOF
OR_test(
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI<0.9))) 
)

# XL MISSENSE
OR_test(
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$mis_z>=3.09))),
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$mis_z<3.09))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z>=3.09))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z<3.09)))
)

# XL LOF
OR_test(
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$pLI>=0.9))),
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$pLI<0.9))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI>=0.9))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI<0.9)))
)


# TRYING WITH OMIM NON-LETHAL VS LETHAL
# AR MISSENSE
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$mis_z<3.09)))
)

# AR LOF
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AR"&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AR"&universe_df$pLI<0.9)))
)

# AR,AD MISSENSE
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$mis_z<3.09)))
)

# AR,AD LOF
OR_test(
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AR,AD"&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AR,AD"&universe_df$pLI<0.9)))  
)

# AD MISSENSE
OR_test(
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$mis_z>=3.09))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$mis_z<3.09))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09)))
)

# AD LOF
OR_test(
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which((universe_df$Inheritance_pattern=="AD"&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI>=0.9))),
  length(which((universe_df$lethal_inheritance=="AD"&universe_df$pLI<0.9))) 
)

# XL MISSENSE
OR_test(
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$mis_z>=3.09))),
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$mis_z<3.09))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z>=3.09))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$mis_z<3.09)))
)

# XL LOF
OR_test(
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$pLI>=0.9))),
  length(which(((universe_df$Inheritance_pattern=="XLr"|universe_df$Inheritance_pattern=="XLd"|universe_df$Inheritance_pattern=="XLd,XLr")&universe_df$human_lethal_B=="N"&universe_df$pLI<0.9))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI>=0.9))),
  length(which(((universe_df$lethal_inheritance=="XLr"|universe_df$lethal_inheritance=="XLd"|
                   universe_df$lethal_inheritance=="XLd,XLr"|universe_df$lethal_inheritance=="XL")&universe_df$pLI<0.9)))
)
