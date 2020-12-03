g.rates <- function(mouse,growth.data){
  df <- growth.data[growth.data$mouse==mouse,]
  ratios <- c(df$vol,NA)/c(NA,df$vol)
  t1 <- c(NA,1,3,6,9,13,NA)
  t2 <- c(NA,3,6,9,13,15,NA)
  #times <- c(NA,"1-3","3-6","6-9","9-13","13-15",NA)
  t1 <- t1[!is.na(ratios)]
  t2 <- t2[!is.na(ratios)]
  ratios <- ratios[!is.na(ratios)]
  
  data.frame(ratios=ratios,t1=t1, t2=t2, mouse=rep(mouse,length(t1)),start=rep(df$start[1],length(t1)),treatment = rep(df$treatment[1],length(t1)))
  
}


load_data  <- function(path="/home/richard/Documents/03_CTLsGoesNative/01_draft_0/01_processed_data/00_formatted/"){
  p1 <- paste(path,"formattedGrowth.csv",sep="")
  p2 <- paste(path,"formattedIntravital.csv",sep="")
  d1 <- read.csv(p1)
  d2 <- read.csv(p2)
 list(gr=d1,iv=d2)
}

get.g <- function(treatment=NA,path="/home/richard/Documents/03_CTLsGoesNative/01_draft_0/01_processed_data/00_formatted/formattedGrowth.csv",ignore_d7=TRUE,post_CTLs=0){
  growth.data <- read.csv(path)
  if(ignore_d7==TRUE) growth.data <- growth.data[!(growth.data$start=="d7.s"),]
  if(post_CTLs==1){
    bool3 <- (growth.data$start=="d3.s") & (growth.data$day<3)
    bool7 <-   (growth.data$start=="d7.s") & (growth.data$day<7)
    growth.data$vol[bool3|bool7] <- NA
  }
  if(post_CTLs==2){
    bool3 <- (growth.data$start=="d3.s") & (growth.data$day>3)
    bool7 <-   (growth.data$start=="d7.s") & (growth.data$day>7)
    growth.data$vol[bool3|bool7] <- NA
  }
  mice <- unique(growth.data$mouse)
  mdf <- do.call(rbind,lapply(mice,g.rates,growth.data))
  mdf$g <- log(mdf$ratios)/(mdf$t2-mdf$t1)
  if(!is.na(treatment)) mdf <- mdf[mdf$treatment==treatment,]
  mdf
}
gr <- get.g(ignore_d7 = FALSE,post_CTLs = TRUE)
get.iv <- function(treatment=NA,path="/home/richard/Documents/03_CTLsGoesNative/01_draft_0/01_processed_data/00_formatted/formattedIntravital.csv"){
  vwin <- 0.35^2*0.1
  rawdf <- read.csv(path,stringsAsFactors = F)
  scale.df <- data.frame(TC=rawdf$TC, CTL=rawdf$CTL,Tden=rawdf$TC/vwin, ET=rawdf$CTL/rawdf$TC, killing=(rawdf$uTC_A+rawdf$cTC_A)/rawdf$CTL, TCa = (rawdf$uTC_A+rawdf$cTC_A)/rawdf$TC, CTLm=rawdf$CTL_M/rawdf$CTL,
                         TCm=rawdf$TC_M/rawdf$TC,CTLa = rawdf$CTL_A/rawdf$CTL,day=rawdf$day,ID=rawdf$ID,mouse=rawdf$mouse,pos=rawdf$pos)
  scale.df[,c("killing","TCa","CTLm","TCm","CTLa")] <- 60*24*scale.df[,c("killing","TCa","CTLm","TCm","CTLa")]/rawdf$min
  if(!is.na(treatment)) scale.df <- scale.df[scale.df$ID==treatment,]
  scale.df
}

