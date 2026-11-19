path <- normalizePath('cause_specific/library/quantreg',mustWork=TRUE)
files <- c(file.path(path,'DESCRIPTION'),list.files(file.path(path,'R'),full.names=TRUE),
           list.files(file.path(path,'libs'),recursive=TRUE,full.names=TRUE))
saveRDS(list(repair='quantreg-6.1-crqfnb-n1-initialisation',package_path=path,
  hashes=tools::md5sum(files)), 'cause_specific/library/repair_manifest.rds')
