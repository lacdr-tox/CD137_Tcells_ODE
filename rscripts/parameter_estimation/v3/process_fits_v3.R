
source("model_functions_v3.R")
require(parallel)


evaluate_fits <- function(ff,fit_dir,d7s=FALSE){ 
  fit <- readRDS(paste(fit_dir,ff,sep="/"))
  x <- unlist(strsplit(ff,split="_"))
  treatment <- x[1]
  id <- unlist(strsplit(x[length(x)],split=".rds"))[1]
  PAR <- fit$optim$bestmem
  fn.control <- fit$fn.control
 
  
  par <- fn.control$unscale(PAR,fn.control$fitted_pars)
  par <- c(par,c(treatment=treatment,id=id))
  out <- run.mod(PAR,fn.control,d7s=d7s)
  df <- iv.sim(id,treatment,PAR,fn.control,out)
  out$g <- par["g"]
  df$g <- par["g"]
  
  out$QSS <- par["QSS"]
  out$s1 <- as.numeric(par["s1"])
  out$id <- as.numeric(id)
  out$treatment <- treatment
  df$QSS <- par["QSS"]
  df$d7s=d7s
  out$d7s=d7s
  
  errordf <- data.frame(treatment=treatment,error=as.numeric(fit$optim$bestval),d7s=d7s,QSS=par["QSS"],g=par["g"],id=id)
  
  par <- par[order(names(par))]
  x <- list(out=out,df=df,par=par,errordf =errordf)
  return(x)
}

## load a saved evolutionary algorithm result, evaulate error for all members, return parameters and error. 
evaluate_members <- function(ff,fit_dir){
  fit <- readRDS(paste(fit_dir,ff,sep="/"))
  x <- unlist(strsplit(ff,split="_"))
  treatment <- x[1]
  id <- unlist(strsplit(x[length(x)],split=".rds"))[1]
  pop <- fit$member$pop
  fn.control <- fit$fn.control
  
  ## find out value of g used 
  PAR <- fit$optim$bestmem
  par <- fn.control$unscale(PAR,fn.control$fitted_pars)
  g <- par["g"]
    
  errors <- unlist(mclapply(1:nrow(pop), function(i){
    PAR <- pop[i,]
    err <- opt.t(PAR,fn.control)
    return(err)},mc.cores = 12))


  x <- list(pop=pop,errors=errors, treatment=treatment,id=id,g=g)
  return(x)
}

## get error for all members of a population
evaluate_population <- function(pop,fn.control){
  errors <- sapply(1:nrow(pop), function(i){
    PAR <- pop[i,]
    err <- opt.t(PAR,fn.control)
    return(err)})

 
  return(errors)
}

process.population <- function(pop){
  x <- data.frame(pop$pop)
  x$errors <- pop$errors
  x$g <- pop$g
  x$id <- pop$id
  x$treatment <- pop$treatment
  return(x)
}

