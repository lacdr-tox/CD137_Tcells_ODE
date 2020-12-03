## usage > get_volgrowth(source="/home/richard/Documents/03_CTLsGoesNative/01_draft_2/01_processed_data/00_formatted/formattedGrowthAll.csv")
g.rates.2 <- function(mouse,growth.data){
  df <- growth.data[growth.data$mouse==mouse,]
  ratios <- df$vol[-1]/df$vol[-length(df$vol)] ## tumour volume ratio between timepoints
  t1 <- df$day[-nrow(df)] # intervals beginning...
  t2 <- df$day[-1] # ... ending
  tmid <- (t1+t2)/2 ## midpoint time of measurements
  v1 <- df$vol[-nrow(df)] # intervals beginning...
  v2 <- df$vol[-1] # ... ending
  vmid <- (v1+v2)/2 ## midpoint time of measurements
  
  dt <- t2-t1 ## lenght of intervals
  g <- log(ratios)/dt ## growth rate estimate
  data.frame(ratios=ratios, g=g, dt=dt, t1=t1, tmid=tmid,t2=t2,v1=v1,v2=v2,vmid=vmid, mouse=rep(mouse,length(t1)),start=rep(df$start[1],length(t1)),treatment = rep(df$treatment[1],length(t1)))
}

get_volgrowth <- function(source, clean=TRUE){ ## there should be a couple of NaNs in the data where mouse 50 was i suppose unmeasured. 1 Inf corresponding to "0" tumour vol on day 1, mouse 56
  growth.data <- read.csv(source)
  mice <- unique(growth.data$mouse)
  df <- do.call(rbind,lapply(mice,g.rates.2,growth.data))
  if(clean==TRUE) df <- df[is.finite(df$g),]
  return(df)
}

