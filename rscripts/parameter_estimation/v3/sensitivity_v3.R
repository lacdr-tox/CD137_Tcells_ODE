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