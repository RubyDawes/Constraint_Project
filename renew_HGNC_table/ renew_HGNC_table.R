#run at beginning of project to update HGNChelper hgnc.table with new gene names- saved in output/Data/hgnc.table.rda so no need to re-run

source(system.file("hgncLookup.R", package = "HGNChelper"))
##You should save this if you are going to use it multiple times,
##then load it from file rather than burdening HGNC's servers.
save(hgnc.table, file="output/Data/hgnc.table.rda", compress="bzip2")