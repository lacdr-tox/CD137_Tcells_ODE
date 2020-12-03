setwd("~/Documents/03_CTLsGoesNative/01_draft_2/01_processed_data/06_image_analysis_full/d7s_mitosis/")

process_file <- function(f){
  df <- read.csv(f)
  if(nrow(df)==0) df <- data.frame(X=NA,Y=NA,Slice=NA,Frame=NA)
  id <- strsplit(f,split="_")[[1]][1:3]
  df <- df[,names(df)%in%c("X","Y","Slice","Frame")]
  df$mouse <- id[1]
  df$day <- id[2]
  df$pos <- strsplit(id[3],split=".csv")[[1]][1]
  df$total <- nrow(df[!is.na(df$X),])
  df
}
## measure mitosis
ff <- list.files()
data <- lapply(ff,process_file)
data <- do.call(rbind,data)
data$treatment <- "ACT-only"
data$treatment[data$mouse%in%c("m48","m54","m56")]<-"ACT+mAb"
data$type <- "mitosis"

mitdf <- data

setwd("~/Documents/03_CTLsGoesNative/01_draft_2/01_processed_data/06_image_analysis_full/d7s_CTLs/")
ff <- list.files()
df <- lapply(ff,process_file)
df <- do.call(rbind,df)
df$treatment <- "ACT-only"
df$treatment[df$mouse%in%c("m48","m54","m56")]<-"ACT+mAb"
df$type <- "CTL"

ctldf <- df

## aggregate mitosis and apoptosis
df <- rbind(data,df)
df <- df[order(df$treatment),]
df$Slice <- df$Slice*7
df$X <- df$X*0.7
df$Y <- df$Y*0.7

## quantify tcs
setwd("/home/richard/Documents/03_CTLsGoesNative/cd137_antiproliferative_ctls/data/results/autocounts_TCs/d7s")
get_TCs <- function(f){
  df <- read.csv(f)
  total <- nrow(df)
  id <- strsplit(f,split="_")[[1]][1:3]
  df <- data.frame(mouse = id[1],
                   day = id[2],
                   pos = strsplit(id[3],split=".csv")[[1]][1],
                   TCs=total)
  
  
  df
}
ff <- list.files("~/Documents/03_CTLsGoesNative/01_draft_2/01_processed_data/06_image_analysis_full/d7s_mitosis/")
tcdf <- lapply(ff,get_TCs)
tcdf <- do.call(rbind,tcdf)
tcdf <- tcdf[order(tcdf$mouse),]
tcdf <- tcdf[order(tcdf$pos),]
## quite hacky:
summary <- aggregate(df[c("total")],by=df[c("treatment","mouse","pos","type")],function(x) x[1])
summary <- summary[order(summary$mouse),]
summary <- summary[order(summary$pos),]
s1 <- summary[summary$type=="CTL",]
s2 <- summary[summary$type=="mitosis",]
summary <- s1
summary$mitosis <- s2$total
names(summary)[names(summary)=="total"] <- "CTL"
summary$TCs <- tcdf$TCs
summary <- summary[,c(1,2,3,5,6,7)]
write.csv(summary,"/home/richard/Documents/03_CTLsGoesNative/cd137_antiproliferative_ctls/data/results/d7s_summary.csv")
## endhack

p <- ggplot(summary,aes(x=CTL,y=mitosis*24/TCs,color=treatment,shape=treatment))+
  geom_point()+
  theme_classic()
p

require(ggplot2)
SAVE <- FALSE
df <- df[!is.na(df$X),]
df <- df[df$Frame<30&df$Slice<(16*7),]
p <- ggplot(df,aes(x=X,y=Slice,colour=type))+
  facet_grid(cols=vars(interaction(treatment,mouse)),rows=vars(pos))+
  geom_point(alpha=0.5)+
  scale_colour_viridis_d()+
  theme_classic()+
  scale_x_continuous(expression(X~(paste(mu,m))))+
  scale_y_continuous(expression(Z~(paste(mu,m))))+
  coord_fixed(ratio=1)

if(SAVE==TRUE) ggsave(filename = "/home/richard/Documents/03_CTLsGoesNative/cd137_antiproliferative_ctls/data/results/manual_counts/XZ.jpg",height=8,width = 8,units="in")

p <- ggplot(df,aes(x=Slice,y=Y,colour=type))+
  facet_grid(cols=vars(interaction(treatment,mouse)),rows=vars(pos))+
  geom_point(alpha=0.5,show.legend = F)+
  scale_colour_viridis_d()+
  theme_classic()+
  scale_x_continuous(expression(Z~(paste(mu,m))))+
  scale_y_continuous(expression(Y~(paste(mu,m))))+
  coord_fixed(ratio=1)


if(SAVE==TRUE)  ggsave(filename = "/home/richard/Documents/03_CTLsGoesNative/cd137_antiproliferative_ctls/data/results/manual_counts/ZY.jpg",height=8,width = 8,units="in")
p <- ggplot(df,aes(x=X,y=Y,colour=type))+
  facet_grid(cols=vars(interaction(treatment,mouse)),rows=vars(pos))+
  geom_point(alpha=0.5,show.legend = F)+
  scale_colour_viridis_d()+
  theme_classic()+
  scale_x_continuous(expression(X~(paste(mu,m))))+
  scale_y_continuous(expression(Y~(paste(mu,m))))+
  coord_fixed(ratio=1)

if(SAVE==TRUE)  ggsave(filename = "/home/richard/Documents/03_CTLsGoesNative/cd137_antiproliferative_ctls/data/results/manual_counts/XY.jpg",height=8,width = 8,units="in")