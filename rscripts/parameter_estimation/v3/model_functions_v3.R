require(deSolve)

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

run.mod <- function(PAR,fn.control,id=1,tmax=15,d7s=FALSE){
  times <- seq(0,tmax,0.1)
  fit.pars <- fn.control$unscale(PAR,fn.control$fitted_pars)
  #print(fit.pars)
  if(d7s==TRUE) fit.pars["t_in"] <- 7
  state <- c(E=0,I=0,R=0,Tc=as.numeric(fit.pars['T0']),Tq=0)
  out <- ode(state,times,tumour.dynamics, fit.pars)
  out <- as.data.frame(out)
  out$id <- id
  return(out)
}

## run model with supplied settings. return all model output
run.mod.sensitivity <- function(PAR,fn.control,id=1,tmax=15,d7s=FALSE,vary="ke",ratio=2){
  times <- seq(0,tmax,0.1)
  fit.pars <- fn.control$unscale(PAR,fn.control$fitted_pars)
  #print(fit.pars)
  if(d7s==TRUE) fit.pars["t_in"] <- 7
  fit.pars[vary] <- fit.pars[vary]*ratio
  state <- c(E=0,I=0,R=0,Tc=as.numeric(fit.pars['T0']),Tq=0)
  out <- ode(state,times,tumour.dynamics, fit.pars)
  out <- as.data.frame(out)
  out$id <- id
  out$ratio<-ratio
  out$vary <- vary
  return(out)
}






iv.sim <- function(id,treatment,PAR,fn.control,sim){
  fit.pars <- fn.control$unscale(PAR,fn.control$fitted_pars)
  sim <- as.data.frame(sim)
  df <- data.frame(time=sim["time"],
                   ET=sim["E"]/rowSums(sim[c("Tc","Tq")]),
                   TCm=fit.pars["g"]*sim["Tc"]/rowSums(sim[c("Tc","Tq")]),
                   killing=fit.pars["ke"],
                   CTLm=sim["I"],
                   CTLa=sim["R"])
  names(df) <- c("time","ET","TCm","killing","CTLm","CTLa")
  if(fit.pars["QSS"]==1){
    I <- fit.pars["ki"]*sim["E"]/rowSums(sim[c("Tc","Tq")])
    df$CTLm<-unlist(I)
  }
  df$id<-id
  df$treatment<- treatment
  return(df)
}

if(FALSE){## example of use
  fn.control <- list(mdf=NA,idf=NA,unscale = function(fit.pars){
    scalars <- c(5000,50,50,50,50,50)
    fit.pars <- fit.pars*scalars
    names(fit.pars) <- c("s1","ki","kr","ke","kq","dq")
    default_pars<-c(const_inf=TRUE,QSS=TRUE,t_in=3,T0=1200,s1 = 50, ki=  0, kr=0, dr = 0,di = 0, g=0.7,  ke=1/3, pow=2/3, kq=10,  dq=2,  t50=8000)
    pars <- c(fit.pars,default_pars[!names(default_pars)%in%names(fit.pars)])
    return(pars)
  })
  PAR<- c(0.131113960, 0.080232033, 0.053239659, 0.012532663, 0.738659715, 0.009738078)
  out <- run.mod(PAR,fn.control)
  df <- iv.sim(1,"agv",PAR,fn.control,out)
  require(ggplot2)
  out <- reshape2::melt(out,id.vars=c("time","id"))
  out <- out[is.finite(log(out$value)),]
  p <- ggplot(out,aes(x=time,y=value))+
    facet_grid(rows=vars(variable),scales="free")+
    geom_line()+
    scale_y_log10()
  plot(p)
}



