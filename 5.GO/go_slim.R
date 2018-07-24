library(GSEABase)

myIds <- c("GO:0016564", "GO:0003677", "GO:0004345", "GO:0008265",
           "GO:0003841", "GO:0030151", "GO:0006355", "GO:0009664",
           "GO:0006412", "GO:0015979", "GO:0006457", "GO:0005618",
           "GO:0005622", "GO:0005840", "GO:0015935", "GO:0000311")
myCollection <- GOCollection(myIds)
fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
slim <- getOBOCollection(fl)
goSlim(myCollection, slim, "MF")
