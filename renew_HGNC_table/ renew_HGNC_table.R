source(system.file("hgncLookup.R", package = "HGNChelper"))
##You should save this if you are going to use it multiple times,
##then load it from file rather than burdening HGNC's servers.
save(hgnc.table, file="hgnc.table.rda", compress="bzip2")