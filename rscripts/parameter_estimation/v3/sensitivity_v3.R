vary_par <- function(ff,fit_dir,d7s=FALSE,vary="ke",range=c(0,1,3,4)){
  fit <- readRDS(paste(fit_dir,ff,sep="/"))
  x <- unlist(strsplit(ff,split="_"))
  treatment <- x[1]
  id <- unlist(strsplit(x[length(x)],split=".rds"))[1]
  PAR <- fit$optim$bestmem
  fn.control <- fit$fn.control
  
  
  par <- fn.control$unscale(PAR,fn.control$fitted_pars)
  par <- c(par,c(treatment=treatment,id=id))
  out <- do.call(rbind,lapply(range, function(r) run.mod.sensitivity(PAR,fn.control,id=id,tmax=30,vary=vary,ratio=r)))
  out$g <- par["g"]
  out$d7s=d7s
  out$treatment <- treatment
  out
}


## run model with supplied settings. return value of objective function
get.error <- function(PAR,fn.control,vary="ke",ratio=2){
  #identify par to vary
  PAR[fn.control$fitted_pars==vary] <- PAR[fn.control$fitted_pars==vary]*ratio
  err <- opt.t(PAR,fn.control)
  return(err)
}

vary_par_error <- function(ff,fit_dir,vary="ke",range=seq(0.9,1.1,0.01)){
  fit <- readRDS(paste(fit_dir,ff,sep="/"))
  x <- unlist(strsplit(ff,split="_"))
  treatment <- x[1]
  PAR <- fit$optim$bestmem
  fn.control <- fit$fn.control
  errors <- sapply(range, function(r) get.error(PAR,fn.control,vary=vary,ratio=r))
  out <- data.frame(errors = errors, range = range, varied = vary, treatment = treatment)
  return(out)
}

sensitivity.analysis <- function(ff, fit_dir,range){
  fit <- readRDS(paste(fit_dir,ff,sep="/"))
  pars <- fit$fn.control$fitted_pars
  
  out <- mclapply(pars, function(pp) vary_par_error(ff,fit_dir,vary=pp,range=range),mc.cores=length(pars))
  out <- do.call(rbind,out)
  return(out)
}
