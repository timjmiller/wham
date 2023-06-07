verify_version = function(model){
  if(!"wham" %in% .packages()) stop("wham library is not loaded")
  wham_commit <- packageDescription("wham")$GithubSHA1
  wham_commit <- ifelse(is.null(wham_commit), "local install", paste0("Github (timjmiller/wham@", wham_commit, ")")) 
  TMB_commit <- packageDescription("TMB")$GithubSHA1
  TMB_commit <- ifelse(is.null(TMB_commit), "local install", paste0("Github (kaskr/adcomp@", TMB_commit, ")")) 
  # print(model$wham_version)
  # print(model$TMB_version)
  # print(model$wham_commit)
  # print(model$TMB_commit)
  # print(wham_commit)
  # print(TMB_commit)
  if(wham_commit !=  model$wham_commit) {
    if(model$wham_commit != "local install") {
      stop(paste0("your wham version: \n ", wham_commit, 
        " \n is not the same as that used to fit your model: \n", 
        model$wham_commit, " \n Load the correct version or refit the model with this version of wham. \n", 
        " If necessary install the right version using \n",
        "devtools::install_github('timjmiller/wham', dependencies=TRUE, ref=", model_commit, ") \n"))
    } else {
      stop(paste0("your wham package is not from Github, but that used to fit your model is: ", model$wham_commit))
    }
  }
  if(TMB_commit !=  model$TMB_commit) {
    if(model$TMB_commit != "local install") {
      stop(paste0("your TMB version: \n ", TMB_commit, 
        " \n is not the same as that used to fit your model: \n", 
        model$TMB_commit, " \n Load the correct version or refit the model with this version of TMB. \n"))
    } else {
      stop(paste0("your TMB package is not from Github, but that used to fit your model is: ", model$TMB_commit))
    }
  }
}
