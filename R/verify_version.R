verify_version = function(model){
  if(!"wham" %in% .packages()) stop("wham library is not loaded")
  if(!"TMB" %in% .packages()) stop("TMB library is not loaded")
  wham_info <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package=="wham") %>% dplyr::select(loadedversion, source) %>% unname
  wham_version <- paste0(wham_info, collapse=" / ")
  model_commit <- strsplit(model$wham_version, "@")[[1]][2]
  model_commit <- substr(model_commit,1,7)
  if(wham_version !=  model$wham_version) {
    stop(paste0("your loaded wham version: \n ", wham_version, 
    " \n is not the same as that used to fit your model: \n", 
    model$wham_version, " \n Load the correct version. If necessary install the right version using \n",
    "devtools::install_github('timjmiller/wham', dependencies=TRUE, ref=", model_commit, ") \n"))
  }
  TMB_info <- sessioninfo::package_info() %>% as.data.frame %>% dplyr::filter(package=="TMB") %>% dplyr::select(loadedversion, source) %>% unname
  TMB_version <- paste0(TMB_info, collapse=" / ")
  model_commit <- strsplit(model$TMB_version, "@")
  model_commit <- substr(model_commit,1,7)
  if(TMB_version !=  model$TMB_version) {
    stop(paste0("your loaded TMB version: \n ", TMB_version, 
    " \n is not the same as that used to fit your model: \n", 
    model$TMB_version, " \n Load the correct version. If necessary install the right version using \n",
    "devtools::install_github('kaskr/adcomp/TMB', dependencies=TRUE, ref=", model_commit, ") \n"))
  }
}
