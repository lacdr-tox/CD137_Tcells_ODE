
source("model_functions_v3.R")



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


