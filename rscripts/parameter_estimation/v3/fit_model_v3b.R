



opt.t <- function(fit.pars,fn.control){
  require(deSolve)
  idf <- fn.control$idf
  mdf <- fn.control$mdf
  idf_d7 <- fn.control$idf_d7
  fitted_pars <- fn.control$fitted_pars
  fit.pars <- fn.control$unscale(fit.pars,fitted_pars)
  
  s1.func <- function(t,t50){
    if(is.na(t50)) return(1)
    else return(1-1/(1+exp(-(t-t50))))
  }
  
  
  tumour.dynamics <- function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      Tt <- Tc+Tq
      dI <- 0
      if(QSS==1) I <- 0
      if(const_inf==1) dE <- as.numeric(t>t_in)*(s1*s1.func(t,(t_in-3+t50))*Tt + (I  - R)*E)
      if(const_inf==0) dE <- as.numeric(t>t_in)*(s1*s1.func(t,(t_in-3+t50))*Tt^(2/3) + (I  - R)*E)
      if(QSS==0) dI <- ki*E/Tt-di*I
      if(QSS==1) dE <- dE+E*ki/Tt
      dR <- kr*E/Tt-dr*R 
      dTc <- g*Tc -(kq+ke)*E*Tc/Tt +dq*Tq
      dTq <- (kq*Tc-ke*Tq)*E/Tt - dq*Tq
      list(c(dE,dI,dR,dTc,dTq)) 
    })
  }
  
  vol.errs <- function(i){
    Tt1 <- sum(out[out$time==mdf$t1[i],c("Tc","Tq")])
    Tt2 <- sum(out[out$time==mdf$t2[i],c("Tc","Tq")])
    dt <- mdf$t2[i]-mdf$t1[i]
    g <- log(Tt2/Tt1)/dt
    err <- (g - mdf$g[i])
    return(err)
  }
  iv.errs <- function(i){
    sim <- out[out$time==idf$day[i],]
    I <- sim["I"]
    if(fit.pars["QSS"]==1){
      I <- fit.pars["ki"]*sim["E"]/rowSums(sim[c("Tc","Tq")])
    }
    meas <- unlist(c(sim["E"]/sum(sim[c("Tc","Tq")]), ##ET RATIO (NOTE THIS IS ONLY FOR ONE TIMEPOINT  (i))
                     fit.pars["g"]*sim["Tc"]/sum(sim[c("Tc","Tq")]), ##GROWTH RATE
                     fit.pars["ke"], ##KILLING RATE
                     I, 
                     sim["R"]))
    names(meas) <- c("ET","TCm","killing","CTLm","CTLa") ## NOTE THAT THESE DO NOT HAVE TO MATCH THE ORDER OF idf, BUT THE NAMES MUST MATCH EXACTLY
    err <- idf[i,names(meas)]-meas
    return(err)
  }
  

  
  state <- c(E=0,I=0,R=0,Tc=as.numeric(fit.pars['T0']),Tq=0)
  times <- seq(0,15,0.1)
  
  mdf_all <- mdf ## we will run the model as usual but just perform the old switcheroo on dataset names to use the funcs above for the d7s data
  mdf <- mdf[mdf$start=="d3.s",]
  err <- tryCatch({
    out <- ode(state,times,tumour.dynamics, fit.pars)
    if(nrow(out)!=length(times)) thiswillleadtoerrors
    out <- data.frame(out)
    v.e.3 <- sapply(1:nrow(mdf), vol.errs)
    i.e <- lapply(1:nrow(idf),iv.errs)
    
    mdf <- mdf_all[mdf_all$start=="d7.s",]
    fit.pars["t_in"]<-7
    out <- ode(state,times,tumour.dynamics, fit.pars)
    if(nrow(out)!=length(times)) thiswillleadtoerrors
    out <- data.frame(out)
    v.e.7 <- sapply(1:nrow(mdf), vol.errs)
    
    ## now include our own manually counted data
    sim <- out[out$time==10,]
    sET <- as.numeric(sim["E"]/sum(sim[c("Tc","Tq")])) ##ET RATIO 
    errET <- idf_d7$ET-sET
    sMIT <- as.numeric(fit.pars["g"]*sim["Tc"]/sum(sim[c("Tc","Tq")])) ##GROWTH RATE
    errMIT <- idf_d7$TCm-sMIT

    errs <- as.numeric(c(v.e.3,v.e.7,unlist(i.e),errET,errMIT))
    n.errs <- length(errs)
    n.exp <- nrow(mdf_all)+5*nrow(idf)+2*nrow(idf_d7)
    n.missing <- n.exp-n.errs
    err <- (sqrt(mean(errs^2,na.rm = T)) + 10^5*n.missing)
    err
    },
    error = function(cond) {return(424242424242)})
  
  return(err)
}



g.rates <- function(mouse,growth.data){
  df <- growth.data[growth.data$mouse==mouse,]
  ratios <- c(df$vol,NA)/c(NA,df$vol)
  t1 <- c(NA,1,3,6,9,13,NA)
  t2 <- c(NA,3,6,9,13,15,NA)
  #times <- c(NA,"1-3","3-6","6-9","9-13","13-15",NA)
  t1 <- t1[!is.na(ratios)&is.finite(ratios)]
  t2 <- t2[!is.na(ratios)&is.finite(ratios)]
  ratios <- ratios[!is.na(ratios)&is.finite(ratios)]
  start <- rep(df$start[1],length(ratios))
  data.frame(ratios=ratios,t1=t1, t2=t2, mouse=rep(mouse,length(t1)),treatment = rep(df$treatment[1],length(t1)),start=start)
  
}

get.g <- function(treatment=NA,path="mouse_data/formattedGrowth.csv"){
  growth.data <- read.csv(path)
  #growth.data <- growth.data[!(growth.data$start=="d7.s"),]
  mice <- unique(growth.data$mouse)
  mdf <- do.call(rbind,lapply(mice,g.rates,growth.data))
  mdf$g <- log(mdf$ratios)/(mdf$t2-mdf$t1)
  if(!is.na(treatment)) mdf <- mdf[mdf$treatment==treatment,]
  mdf 
}

get.iv <- function(treatment=NA,path="mouse_data/formattedIntravital.csv"){
  rawdf <- read.csv(path,stringsAsFactors = F)
  scale.df <- data.frame(ET=rawdf$CTL/rawdf$TC, killing=(rawdf$uTC_A+rawdf$cTC_A)/rawdf$CTL, TCa = (rawdf$uTC_A+rawdf$cTC_A)/rawdf$TC, CTLm=rawdf$CTL_M/rawdf$CTL,
                         TCm=rawdf$TC_M/rawdf$TC,CTLa = rawdf$CTL_A/rawdf$CTL,day=rawdf$day,ID=rawdf$ID,mouse=rawdf$mouse,pos=rawdf$pos)
  scale.df[,c("killing","TCa","CTLm","TCm","CTLa")] <- 60*24*scale.df[,c("killing","TCa","CTLm","TCm","CTLa")]/rawdf$min
  if(!is.na(treatment)) scale.df <- scale.df[scale.df$ID==treatment,]
  scale.df
}

get.d7.data <- function(treatment=NA,path="mouse_data/d7s_summary.csv"){
  rawdf <- read.csv(path,stringsAsFactors = F)
  df <- data.frame(treatment=rawdf$treatment, ET=rawdf$CTL/rawdf$TCs, TCm=24*rawdf$mitosis/rawdf$TCs) 
  if(treatment=="mAB") return(df[df$treatment=="ACT+mAb",])
  if(treatment=="control") return(df[df$treatment=="ACT-only",])
}

fit_model <- function(treatment,CI=FALSE,QSS=FALSE,set_growth=0.7,fitted_pars=c("s1","ki","kr","ke","kq","dq"),save_dir="test_output",save_name="A",gens=200,pop=100,ncore=8){
    ## build the function control, can modify the values (BUT NOT THE NAMES) of default_pars
    ## scalars must match the par.names and is used to set ranges on the fitted pars.
    ## evolutionary algorithm benefits if the values of the fitted parameters have roughly the same size, which is why we do this scaling
    ## length of upper and lowe must also be matched manually
    ## SET QSS TRUE/FALSE for QUASISteadyState, CI for constant infiltration of CTLs into tumour, CI+false for variable rate per volume
    mdf = get.g(treatment)
    idf=get.iv(treatment)## to get around the dependence issue, we could here select only 1 position per mouse
    idf_d7=get.d7.data(treatment)
    fn.control <- list(mdf=mdf,idf=idf,idf_d7=idf_d7,fitted_pars=fitted_pars,unscale = function(fit.pars,fitted_pars){
      scalars <- rep(50,length(fitted_pars))
      scalars[fitted_pars=="s1"] <- scalars[fitted_pars=="s1"]*0.1 
      if(CI==TRUE) scalars[fitted_pars=="s1"] <- scalars[fitted_pars=="s1"]*100
      scalars[fitted_pars=="g"] <- scalars[fitted_pars=="g"]*0.02 
      scalars[fitted_pars=="pow"] <- scalars[fitted_pars=="pow"]*0.02 
      fit.pars <- fit.pars*scalars
      names(fit.pars) <- fitted_pars
      default_pars<-c(const_inf=CI,QSS=QSS,t_in=3,T0=1200,s1 = 50, ki=  0, kr=0, dr = 0,di = 0, g=set_growth,  ke=1/3, pow=2/3,  kq=10,  dq=2,  t50=8000)
      pars <- c(fit.pars,default_pars[!names(default_pars)%in%names(fit.pars)])
      return(pars)
    })
    upper <- rep(1,length(fitted_pars))
    lower <- rep(0,length(fitted_pars))
    require(DEoptim)
  
    #opt <- DEoptim(opt.t,lower=lower,upper=upper, control=list(NP=100,itermax=100),fn.control=fn.control)
    source("deopt2.R")
    opt <- DEoptim2(opt.t,lower=lower,upper=upper, control=list(NP=pop,itermax=gens,parallelType=1,limitCores=ncore),fn.control=fn.control)
    opt$fn.control <- fn.control
    require(stringr)
    fit_id <- length(list.files(save_dir))
    fname <- paste(save_dir,"/",treatment,"_",save_name,"_",str_pad(fit_id,width = 3,pad = 0),".rds",sep="")
    print(fn.control$unscale(opt$optim$bestmem,fitted_pars))
    saveRDS(opt,fname)
}




