# False trends in Landsat-based assessment of vegetation greening in world's cold regions
# due to increasing images availability over time

# Arthur Bayle, Simon Gascoin, Logan T. Berner and Philippe Choler

library(terra)
library(tidyr)
library(dplyr)
library(viridis)
library(tibble)
library(ggplot2)
library(data.table)
library(mblm)
library(randomForest)
library(caret)
setwd("D:/")
'%!in%' <- function(x,y)!('%in%'(x,y))

# ---- I. FIGURE 1 ----
# a. Prepare data ----
# Load DEM
DEM = rast("BIS/DATA/DEM/DEM_EUALPS.tif")

# Number of observations
OBS = list.files("BIS2/DATA/OBS/",pattern="NumberObs",full.names=T)
OBSr = rast(OBS)

SPL_OUT = vect("BIS2/DATA/SAMPLES/ALL_OUT.shp")
SPL_OVL = vect("BIS2/DATA/SAMPLES/ALL_OVL.shp")

# Extract elevation on samples
DEM_OUT = extract(DEM,SPL_OUT)
DEM_OVL = extract(DEM,SPL_OVL)

# Extract observations on samples
OBS_OUT = extract(OBSr,SPL_OUT)
OBS_OVL = extract(OBSr,SPL_OVL)

CSV_OBS_OUT = data.frame(ID = DEM_OUT$ID, DEM = round(DEM_OUT$DEM_EU))
CSV_OBS_OUT = cbind(CSV_OBS_OUT,round(OBS_OUT[,-1]))
colnames(CSV_OBS_OUT)[3:40] = c(1984:2021)
write.table(CSV_OBS_OUT,"BIS2/DATA/CSV/Observations_Outside.csv",quote=F,sep=";",row.names=F)

CSV_OBS_OVL = data.frame(ID = DEM_OVL$ID, DEM = round(DEM_OVL$DEM_EU))
CSV_OBS_OVL = cbind(CSV_OBS_OVL,round(OBS_OVL[,-1]))
colnames(CSV_OBS_OVL)[3:40] = c(1984:2021)
write.table(CSV_OBS_OVL,"BIS2/DATA/CSV/Observations_Overlap.csv",quote=F,sep=";",row.names=F)

# b. Load data ----

OBS_OUT = read.csv("BIS2/DATA/CSV/Observations_Outside.csv",sep=";")
OBS_OVL = read.csv("BIS2/DATA/CSV/Observations_Overlap.csv",sep=";")

APP_OUT = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_OUT[,-c(1,2)], 2, median)))
APP_OUTsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_OUT[,-c(1,2)], 2, sd)))

APP_OVL = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_OVL[,-c(1,2)], 2, median)))
APP_OVLsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_OVL[,-c(1,2)], 2, sd)))

# c. Panel B ----
pdf("BIS4/FIGURES/FIGURE1/PANELS/Figure1_panelA.pdf",width = 7,height = 7)

plot(1,type="n",xlim=c(1984,2021),ylim=c(0,13),yaxs="i",ylab="Useable observations",xlab="Year",xaxs="i")
abline(h=c(2.5,5,7.5,10,12.5,15),lty=3,col="grey90")
abline(v=seq(1985,2020,5),lty=3,col="grey90")

lines(APP_OUT,type="l",lwd=2,col="#7b3294")
polygon(x=c(APP_OUT$YEAR,rev(APP_OUT$YEAR)),y=c(APP_OUT$OBS-APP_OUTsd$OBS,rev(APP_OUT$OBS+APP_OUTsd$OBS)),
        col=adjustcolor("#c2a5cf",alpha.f = 0.1),border = NA)
text(2010,1,"Outside overlapping tiles",cex=1,col="#7b3294")

lines(APP_OVL,type="l",lwd=2,col="#008837",ylim=c(0,17),yaxs="i",ylab="Landsat observations",xlab="Year",xaxs="i")
polygon(x=c(APP_OVL$YEAR,rev(APP_OVL$YEAR)),y=c(APP_OVL$OBS-APP_OVLsd$OBS,rev(APP_OVL$OBS+APP_OVLsd$OBS)),
        col=adjustcolor("#a6dba0",alpha.f = 0.1),border = NA)
text(1992,8,"Within overlapping tiles",cex=1,col="#008837")


arrows(x0 = 1984.2,y0 = 9.5,x1 = 2011,y1 = 9.5,angle = 90,length = 0.1,col="red",lwd=2)
arrows(x1 = 1984.2,y0 = 9.5,x0 = 2011,y1 = 9.5,angle = 90,length = 0.1,col="red",lwd=2)
text(x=1998,y=10,"Landsat 5 TM",col="red")

arrows(x1 = 1999,y0 = 11,x0 = 2021,y1 = 11,angle = 90,length = 0.1,col="black",lwd=2)
for(i in 2003:2020){segments(x0 = i,x1=i+1,y0 = 10.8,y1=11.2,lwd=2)}
text(x=2010,y=10.5,"Landsat 7 ETM+",col="black")

arrows(x1 = 2013,y0 = 12.5,x0 = 2021,y1 = 12.5,angle = 90,length = 0.1,col="orange",lwd=2)
text(x=2017.2,y=12,"Landsat 8 OLI",col="orange")

for(i in 1985:1987){segments(x0 = i,x1=i+1,y0 = 11.8,y1=12.2,lwd=2)}
segments(1985,12,1988,12,lwd=2)
text(x=1990.7,y=12,"SLC-off")

dev.off()


# d. Prepare data ----
OBS_OUT = read.csv("BIS2/DATA/CSV/Observations_Outside.csv",sep=";")
OBS_OVL = read.csv("BIS2/DATA/CSV/Observations_Overlap.csv",sep=";")

DEMdf_OUT = data.frame(DEM = OBS_OUT$DEM,
                       OBSP1 = apply(OBS_OUT[,3:17],1,mean),
                       OBSP2 = apply(OBS_OUT[,18:31],1,mean),
                       OBSP3 = apply(OBS_OUT[,32:40],1,mean))

DEMdf_OUT$CLASS = NA
DEMdf_OUT$CLASS[which(DEMdf_OUT$DEM < 1000)] = "[< 1000m]"
DEMdf_OUT$CLASS[which(DEMdf_OUT$DEM > 1000 & DEMdf_OUT$DEM < 1500)] = "[1000-1500m]"
DEMdf_OUT$CLASS[which(DEMdf_OUT$DEM > 1500 & DEMdf_OUT$DEM < 2000)] = "[1500-2000m]"
DEMdf_OUT$CLASS[which(DEMdf_OUT$DEM > 2000 & DEMdf_OUT$DEM < 2500)] = "[2000-2500m]"
DEMdf_OUT$CLASS[which(DEMdf_OUT$DEM > 2500 & DEMdf_OUT$DEM < 3000)] = "[2500-3000m]"
DEMdf_OUT$CLASS[which(DEMdf_OUT$DEM > 3000)] = "[> 3000m]"

DEMdf_OVL = data.frame(DEM = OBS_OVL$DEM,
                       OBSP1 = apply(OBS_OVL[,3:17],1,mean),
                       OBSP2 = apply(OBS_OVL[,18:31],1,mean),
                       OBSP3 = apply(OBS_OVL[,32:40],1,mean))

DEMdf_OVL$CLASS = NA
DEMdf_OVL$CLASS[which(DEMdf_OVL$DEM < 1000)] = "[< 1000m]"
DEMdf_OVL$CLASS[which(DEMdf_OVL$DEM > 1000 & DEMdf_OVL$DEM < 1500)] = "[1000-1500m]"
DEMdf_OVL$CLASS[which(DEMdf_OVL$DEM > 1500 & DEMdf_OVL$DEM < 2000)] = "[1500-2000m]"
DEMdf_OVL$CLASS[which(DEMdf_OVL$DEM > 2000 & DEMdf_OVL$DEM < 2500)] = "[2000-2500m]"
DEMdf_OVL$CLASS[which(DEMdf_OVL$DEM > 2500 & DEMdf_OVL$DEM < 3000)] = "[2500-3000m]"
DEMdf_OVL$CLASS[which(DEMdf_OVL$DEM > 3000)] = "[> 3000m]"


# e. Panel C1 ----
pdf("BIS5/FIGURES/FIGURE1/PANELS/Figure1_panelC1.pdf",width = 5,height = 7)
boxplot( DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[< 1000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[< 1000m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[1000-1500m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[1000-1500m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[1500-2000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[1500-2000m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[2000-2500m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[2000-2500m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[2500-3000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[2500-3000m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[> 3000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[> 3000m]")],
         
         
         names=c(),outline=F,col=rep(adjustcolor(c("#008837","#7b3294"),alpha.f = 0.5)),
         xaxs="i",xlim=c(1,12),ylim=c(0,17),xaxt="n",main="1984 - 1998",pch=".",boxwex=0.7
         ,yaxt="n")
axis(side = 2, at = seq(0,20,2),las=2)

abline(v=seq(2.5,16.5,2))
abline(h=seq(0,20,1),lty=2,col=adjustcolor("grey",alpha.f = 0.7))

boxplot( DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[< 1000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[< 1000m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[1000-1500m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[1000-1500m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[1500-2000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[1500-2000m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[2000-2500m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[2000-2500m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[2500-3000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[2500-3000m]")],
         
         DEMdf_OVL$OBSP1[which(DEMdf_OVL$CLASS == "[> 3000m]")],
         DEMdf_OUT$OBSP1[which(DEMdf_OUT$CLASS == "[> 3000m]")],
         
         
         names=c(),outline=F,col=rep(adjustcolor(c("#008837","#7b3294"),alpha.f = 1)),
         xaxs="i",xlim=c(1,12),ylim=c(0,17),xaxt="n",ylab="USeable observations",add=T,pch=".",boxwex=0.7,yaxt="n")
dev.off()
# c. Panel C2 ----
pdf("BIS5/FIGURES/FIGURE1/PANELS/Figure1_panelC2.pdf",width = 5,height = 7)

boxplot( DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[< 1000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[< 1000m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[1000-1500m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[1000-1500m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[1500-2000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[1500-2000m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[2000-2500m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[2000-2500m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[2500-3000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[2500-3000m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[> 3000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[> 3000m]")],
         
         
         names=c(),outline=F,col=rep(adjustcolor(c("#008837","#7b3294"),alpha.f = 0.5)),
         xaxs="i",xlim=c(1,12),ylim=c(0,17),xaxt="n",main="1999 - 2012",pch=".",boxwex=0.7,yaxt="n")
axis(side = 2, at = seq(0,20,2),las=2)

abline(v=seq(2.5,16.5,2))
abline(h=seq(0,20,1),lty=2,col=adjustcolor("grey",alpha.f = 0.7))

boxplot( DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[< 1000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[< 1000m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[1000-1500m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[1000-1500m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[1500-2000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[1500-2000m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[2000-2500m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[2000-2500m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[2500-3000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[2500-3000m]")],
         
         DEMdf_OVL$OBSP2[which(DEMdf_OVL$CLASS == "[> 3000m]")],
         DEMdf_OUT$OBSP2[which(DEMdf_OUT$CLASS == "[> 3000m]")],
         
         
         names=c(),outline=F,col=rep(adjustcolor(c("#008837","#7b3294"),alpha.f = 1)),
         xaxs="i",xlim=c(1,12),ylim=c(0,17),xaxt="n",ylab="USeable observations",add=T,pch=".",boxwex=0.7,yaxt="n")
dev.off()
# d. Panel C3 ----
pdf("BIS5/FIGURES/FIGURE1/PANELS/Figure1_panelC3.pdf",width = 5,height = 7)

boxplot( DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[< 1000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[< 1000m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[1000-1500m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[1000-1500m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[1500-2000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[1500-2000m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[2000-2500m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[2000-2500m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[2500-3000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[2500-3000m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[> 3000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[> 3000m]")],
         
         
         names=c(),outline=F,col=rep(adjustcolor(c("#008837","#7b3294"),alpha.f = 0.5)),
         xaxs="i",xlim=c(1,12),ylim=c(0,17),xaxt="n",main="2013 - 2021",pch=".",boxwex=0.7,yaxt="n")
axis(side = 2, at = seq(0,20,2),las=2)

abline(v=seq(2.5,16.5,2))
abline(h=seq(0,20,1),lty=2,col=adjustcolor("grey",alpha.f = 0.7))

boxplot( DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[< 1000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[< 1000m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[1000-1500m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[1000-1500m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[1500-2000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[1500-2000m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[2000-2500m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[2000-2500m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[2500-3000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[2500-3000m]")],
         
         DEMdf_OVL$OBSP3[which(DEMdf_OVL$CLASS == "[> 3000m]")],
         DEMdf_OUT$OBSP3[which(DEMdf_OUT$CLASS == "[> 3000m]")],
         
         
         names=c(),outline=F,col=rep(adjustcolor(c("#008837","#7b3294"),alpha.f = 1)),
         xaxs="i",xlim=c(1,12),ylim=c(0,17),xaxt="n",ylab="USeable observations",add=T,pch=".",boxwex=0.7,yaxt="n")
dev.off()


# ---- II. FIGURE 3 ----
# a. Prepare data ----
DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

# Load pheno
PHENO = read.table("BIS/DATA/PHENO/PHENO_CLUST4.csv",header=T,sep=",")
PARAM = read.table("BIS/DATA/PHENO/OPTIM_CLUST4.csv",header=T,sep=",")

# b. Panel A ----
pdf("BIS5/FIGURES/FIGURE3/PANELS/Figure3_pheno.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))

plot(y=DL.f(PAR = PARAM$cluster2,DOY = DOY.seq),x=DOY.seq,ylim=c(0,1),type="l",lwd=3,col="red",
     ylab="NDVI",xlab="Day of year",xaxs="i",xlim=c(60,360),lty=2)
#polygon(x=c(152,243,243,152),y=c(-0.1,-0.1,1.1,1.1),col=adjustcolor("grey",alpha.f = 0.1),border = NA)
lines(PHENO$X,PHENO$cluster2/1000,ylim=c(0,1),type="l",lwd=3,col="red",ylab="NDVI",xlab="Day of year",xaxs="i")
#polygon(x = c(DOY.seq,rev(DOY.seq)),y=c(DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)-0.1,rev(DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)+0.1)),
#        col=adjustcolor("red",alpha.f = 0.2),border = NA)

points(DL.f(PAR = PARAM$cluster3,DOY = DOY.seq),x=DOY.seq,type="l",lwd=3,col="black",lty=2)
lines(PHENO$X,PHENO$cluster3/1000,type="l",lwd=3,col="black")
#polygon(x = c(DOY.seq,rev(DOY.seq)),y=c(DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)-0.1,rev(DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)+0.1)),
#        col=adjustcolor("black",alpha.f = 0.2),border = NA)

points(DL.f(PAR = PARAM$cluster4,DOY = DOY.seq),x=DOY.seq,type="l",lwd=3,col="cadetblue",lty=2)
lines(PHENO$X,PHENO$cluster4/1000,type="l",lwd=3,col="cadetblue")
#polygon(x = c(DOY.seq,rev(DOY.seq)),y=c(DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)-0.1,rev(DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)+0.1)),
#        col=adjustcolor("cadetblue",alpha.f = 0.2),border = NA)

#text(x=70,y=0.98,"Thermic grasslands (GSL ~ 150)",cex=1,adj=0,col="black")
#text(x=70,y=0.93,"Intermediate grasslands (GSL ~ 125)",cex=1,adj=0,col="red")
#text(x=70,y=0.88,"Snowbeds grasslands (GSL ~ 100)",cex=1,adj=0,col="cadetblue")

dev.off()

# c. Prepare data ----
NDVIyts_CL2 = DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)
NDVIyts_CL3 = DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)
NDVIyts_CL4 = DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)

DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

# Sample NDVI under varying revisit intervals and estimate NDVImax using a retrieving function ---
NSAMP       = 10000         # number of sampling
REVISIT.freq = 8   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]

# NOBSmin is the minimum number of dates for a given revisit interval  
NOBSmin  <- min(unlist(lapply(DO.tmp,function(x) length(x)))) 

# initialize the matrix to store NDVImax estimates  
NDVImax1.tmp_CL2   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL3   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL4   <- matrix(ncol=NOBSmin,nrow=NSAMP)

DOYgs.tmp_CL2      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL3      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL4      <- matrix(ncol=NOBSmin,nrow=NSAMP)


for (i in 1:NOBSmin){
  print(i)
  
  tmp_CL2     <- sapply(1:NSAMP, function(x) NDVIyts_CL2[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL2 = matrix(ncol=length(tmp_CL2),data=tmp_CL2)}
  tmp.list_CL2 <- split(tmp_CL2, seq(ncol(tmp_CL2)))
  
  tmp_CL3     <- sapply(1:NSAMP, function(x) NDVIyts_CL3[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL3 = matrix(ncol=length(tmp_CL3),data=tmp_CL3)}
  tmp.list_CL3 <- split(tmp_CL3, seq(ncol(tmp_CL3)))
  
  tmp_CL4     <- sapply(1:NSAMP, function(x) NDVIyts_CL4[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL4 = matrix(ncol=length(tmp_CL4),data=tmp_CL4)}
  tmp.list_CL4 <- split(tmp_CL4, seq(ncol(tmp_CL4)))
  
  
  DOYgs.f            <- function(x) length(x)
  DOYgs.tmp_CL2[,i]    <- as.numeric(unlist(lapply(tmp.list_CL2,DOYgs.f)))
  DOYgs.tmp_CL3[,i]    <- as.numeric(unlist(lapply(tmp.list_CL3,DOYgs.f)))
  DOYgs.tmp_CL4[,i]    <- as.numeric(unlist(lapply(tmp.list_CL4,DOYgs.f)))
  
  NDVImax.f1          <- function(x) max(x)
  NDVImax1.tmp_CL2[,i]  <- as.numeric(unlist(lapply(tmp.list_CL2,NDVImax.f1)))
  NDVImax1.tmp_CL3[,i]  <- as.numeric(unlist(lapply(tmp.list_CL3,NDVImax.f1)))
  NDVImax1.tmp_CL4[,i]  <- as.numeric(unlist(lapply(tmp.list_CL4,NDVImax.f1)))
  
}  

NOBSrange   <- 1:10

MAX_CL2        <- matrix(NA,length(NOBSrange))
MAX_CL3        <- matrix(NA,length(NOBSrange))
MAX_CL4        <- matrix(NA,length(NOBSrange))

SD_CL2        <- matrix(NA,length(NOBSrange))
SD_CL3        <- matrix(NA,length(NOBSrange))
SD_CL4        <- matrix(NA,length(NOBSrange))

for(i in 1:length(NOBSrange)){ # number of observations (-1)
  SEL          <- which(DOYgs.tmp_CL2==NOBSrange[i],arr.ind=T)
  DF           <- cbind(SEL,NDVImax1=NDVImax1.tmp_CL2[SEL])
  
  MAX_CL2[i] <- (mean(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  MAX_CL3[i] <- (mean(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  MAX_CL4[i] <- (mean(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
  SD_CL2[i] <- (sd(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  SD_CL3[i] <- (sd(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  SD_CL4[i] <- (sd(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
}

# d. Panel B ----
CLU = rast("BIS2/DATA/PREPARE SAMPLE/DistributionClusters.tif")
TOKEEP = rast("BIS2/DATA/PREPARE SAMPLE/FinalSample.tif")
CLU[which(TOKEEP)]
DEM = rast("BIS/DATA/DEM/DEM_EUALPS.tif")

# How many pixels per clusters ?
CLU = mask(CLU,TOKEEP)
# Cluster 1 (removed) : 2074179 (13 %)
# Cluster 2 (red) : 5945676 (37 %)
# Cluster 3 (black) : 6711065 (42 %)
# Cluster 4 (blue) : 1035114 (6 %)

# Elevation per clusters
DEM_CL2 = DEM[which(CLU[] == 2)][,1]
DEM_CL3 = DEM[which(CLU[] == 3)][,1]
DEM_CL4 = DEM[which(CLU[] == 4)][,1]

pdf("BIS4/FIGURES/FIGURE3/PANELS/Figure3_clust_elev.pdf",width = 7,height = 6)
hist(DEM_CL3,xlim=c(1400,3000),breaks=50,col=adjustcolor("black",alpha.f = 0.5),main="",
     yaxs="i",freq=F,xlab="Elevation (m)",ylim=c(0,0.0025),xaxs="i")
hist(DEM_CL2,xlim=c(1400,3200),breaks=100,col=adjustcolor("red",alpha.f = 0.5),add=T,freq=F)
hist(DEM_CL4,xlim=c(1400,3200),breaks=100,col=adjustcolor("cadetblue",alpha.f = 0.8),add=T,freq=F)
box()
dev.off()

# e. Prepare data ----
# Load DEM
DEM = rast("BIS/DATA/DEM/DEM_EUALPS.tif")

# Load samples
SPL_OUT = vect("BIS2/DATA/SAMPLES/CLU3_OUT.shp")
SPL_OVL = vect("BIS2/DATA/SAMPLES/CLU3_OVL.shp")

# Load observations
OBS = list.files("BIS2/DATA/OBS/",pattern="NumberObs",full.names=T)
OBSr = rast(OBS)

# Extract elevation on samples
DEM_OUT = extract(DEM,SPL_OUT)
DEM_OVL = extract(DEM,SPL_OVL)

# Extract observations on samples
OBS_OUT = extract(OBSr,SPL_OUT)
OBS_OVL = extract(OBSr,SPL_OVL)

CSV_OBS_OUT = data.frame(ID = DEM_OUT$ID, DEM = round(DEM_OUT$DEM_EU))
CSV_OBS_OUT = cbind(CSV_OBS_OUT,round(OBS_OUT[,-1]))
colnames(CSV_OBS_OUT)[3:40] = c(1984:2021)
write.table(CSV_OBS_OUT,"BIS2/DATA/CSV/Observations_Black_Outside.csv",quote=F,sep=";",row.names=F)

CSV_OBS_OVL = data.frame(ID = DEM_OVL$ID, DEM = round(DEM_OVL$DEM_EU))
CSV_OBS_OVL = cbind(CSV_OBS_OVL,round(OBS_OVL[,-1]))
colnames(CSV_OBS_OVL)[3:40] = c(1984:2021)
write.table(CSV_OBS_OVL,"BIS2/DATA/CSV/Observations_Black_Overlap.csv",quote=F,sep=";",row.names=F)

# Load DEM
DEM = rast("BIS/DATA/DEM/DEM_EUALPS.tif")

# Load samples
SPL_OUT = vect("BIS2/DATA/SAMPLES/CLU2_OUT.shp")
SPL_OVL = vect("BIS2/DATA/SAMPLES/CLU2_OVL.shp")

# Load observations
OBS = list.files("BIS2/DATA/OBS/",pattern="NumberObs",full.names=T)
OBSr = rast(OBS)

# Extract elevation on samples
DEM_OUT = extract(DEM,SPL_OUT)
DEM_OVL = extract(DEM,SPL_OVL)

# Extract observations on samples
OBS_OUT = extract(OBSr,SPL_OUT)
OBS_OVL = extract(OBSr,SPL_OVL)

CSV_OBS_OUT = data.frame(ID = DEM_OUT$ID, DEM = round(DEM_OUT$DEM_EU))
CSV_OBS_OUT = cbind(CSV_OBS_OUT,round(OBS_OUT[,-1]))
colnames(CSV_OBS_OUT)[3:40] = c(1984:2021)
write.table(CSV_OBS_OUT,"BIS2/DATA/CSV/Observations_Red_Outside.csv",quote=F,sep=";",row.names=F)

CSV_OBS_OVL = data.frame(ID = DEM_OVL$ID, DEM = round(DEM_OVL$DEM_EU))
CSV_OBS_OVL = cbind(CSV_OBS_OVL,round(OBS_OVL[,-1]))
colnames(CSV_OBS_OVL)[3:40] = c(1984:2021)
write.table(CSV_OBS_OVL,"BIS2/DATA/CSV/Observations_Red_Overlap.csv",quote=F,sep=";",row.names=F)

# Load DEM
DEM = rast("BIS/DATA/DEM/DEM_EUALPS.tif")

# Load samples
SPL_OUT = vect("BIS2/DATA/SAMPLES/CLU4_OUT.shp")
SPL_OVL = vect("BIS2/DATA/SAMPLES/CLU4_OVL.shp")

# Load observations
OBS = list.files("BIS2/DATA/OBS/",pattern="NumberObs",full.names=T)
OBSr = rast(OBS)

# Extract elevation on samples
DEM_OUT = extract(DEM,SPL_OUT)
DEM_OVL = extract(DEM,SPL_OVL)

# Extract observations on samples
OBS_OUT = extract(OBSr,SPL_OUT)
OBS_OVL = extract(OBSr,SPL_OVL)

CSV_OBS_OUT = data.frame(ID = DEM_OUT$ID, DEM = round(DEM_OUT$DEM_EU))
CSV_OBS_OUT = cbind(CSV_OBS_OUT,round(OBS_OUT[,-1]))
colnames(CSV_OBS_OUT)[3:40] = c(1984:2021)
write.table(CSV_OBS_OUT,"BIS2/DATA/CSV/Observations_Blue_Outside.csv",quote=F,sep=";",row.names=F)

CSV_OBS_OVL = data.frame(ID = DEM_OVL$ID, DEM = round(DEM_OVL$DEM_EU))
CSV_OBS_OVL = cbind(CSV_OBS_OVL,round(OBS_OVL[,-1]))
colnames(CSV_OBS_OVL)[3:40] = c(1984:2021)
write.table(CSV_OBS_OVL,"BIS2/DATA/CSV/Observations_Blue_Overlap.csv",quote=F,sep=";",row.names=F)

# f. Load data ----

# Load black cluster
OBS_BLACK_OUT = read.csv("BIS2/DATA/CSV/Observations_Black_Outside.csv",sep=";")
OBS_BLACK_OVL = read.csv("BIS2/DATA/CSV/Observations_Black_Overlap.csv",sep=";")

APP_BLACK_OUT = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLACK_OUT[,-c(1,2)], 2, mean)))
APP_BLACK_OUTsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLACK_OUT[,-c(1,2)], 2, sd)))

APP_BLACK_OVL = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLACK_OVL[,-c(1,2)], 2, mean)))
APP_BLACK_OVLsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLACK_OVL[,-c(1,2)], 2, sd)))

# Load red cluster
OBS_RED_OUT = read.csv("BIS2/DATA/CSV/Observations_Red_Outside.csv",sep=";")
OBS_RED_OVL = read.csv("BIS2/DATA/CSV/Observations_Red_Overlap.csv",sep=";")

APP_RED_OUT = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_RED_OUT[,-c(1,2)], 2, mean)))
APP_RED_OUTsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_RED_OUT[,-c(1,2)], 2, sd)))

APP_RED_OVL = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_RED_OVL[,-c(1,2)], 2, mean)))
APP_RED_OVLsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_RED_OVL[,-c(1,2)], 2, sd)))

# Load blue cluster
OBS_BLUE_OUT = read.csv("BIS2/DATA/CSV/Observations_Blue_Outside.csv",sep=";")
OBS_BLUE_OVL = read.csv("BIS2/DATA/CSV/Observations_Blue_Overlap.csv",sep=";")

APP_BLUE_OUT = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLUE_OUT[,-c(1,2)], 2, mean)))
APP_BLUE_OUTsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLUE_OUT[,-c(1,2)], 2, sd)))

APP_BLUE_OVL = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLUE_OVL[,-c(1,2)], 2, mean)))
APP_BLUE_OVLsd = data.frame(YEAR = 1984:2021, OBS = as.numeric(apply(OBS_BLUE_OVL[,-c(1,2)], 2, sd)))

# g. Panel C and D ----
pdf("BIS5/FIGURES/FIGURE3/PANELS/Figure3_ovl_trends.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))
plot(1,type="n",xlim=c(1984,2021),ylim=c(0,10),yaxs="i",ylab="Useable observations",xlab="Year",xaxs="i")
abline(h=c(2.5,5,7.5,10,12.5,15),lty=3,col="grey90")
abline(v=seq(1985,2020,5),lty=3,col="grey90")

lines(APP_BLACK_OVL,type="l",lwd=2,col="black")
polygon(x=c(APP_BLACK_OVL$YEAR,rev(APP_BLACK_OVL$YEAR)),y=c(APP_BLACK_OVL$OBS-APP_BLACK_OVLsd$OBS,rev(APP_BLACK_OVL$OBS+APP_BLACK_OVLsd$OBS)),
        col=adjustcolor("black",alpha.f = 0.1),border = NA)

lines(APP_RED_OVL,type="l",lwd=2,col="red")
polygon(x=c(APP_RED_OVL$YEAR,rev(APP_RED_OVL$YEAR)),y=c(APP_RED_OVL$OBS-APP_RED_OVLsd$OBS,rev(APP_RED_OVL$OBS+APP_RED_OVLsd$OBS)),
        col=adjustcolor("red",alpha.f = 0.1),border = NA)

lines(APP_BLUE_OVL,type="l",lwd=2,col="cadetblue")
polygon(x=c(APP_BLUE_OVL$YEAR,rev(APP_BLUE_OVL$YEAR)),y=c(APP_BLUE_OVL$OBS-APP_BLUE_OVLsd$OBS,rev(APP_BLUE_OVL$OBS+APP_BLUE_OVLsd$OBS)),
        col=adjustcolor("cadetblue",alpha.f = 0.1),border = NA)
dev.off()

pdf("BIS5/FIGURES/FIGURE3/PANELS/Figure3_out_trends.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))

plot(1,type="n",xlim=c(1984,2021),ylim=c(0,10),yaxs="i",ylab="Useable observations",xlab="Year",xaxs="i")
abline(h=c(2.5,5,7.5,10,12.5,15),lty=3,col="grey90")
abline(v=seq(1985,2020,5),lty=3,col="grey90")

lines(APP_BLACK_OUT,type="l",lwd=2,col="black")
polygon(x=c(APP_BLACK_OUT$YEAR,rev(APP_BLACK_OUT$YEAR)),y=c(APP_BLACK_OUT$OBS-APP_BLACK_OUTsd$OBS,rev(APP_BLACK_OUT$OBS+APP_BLACK_OUTsd$OBS)),
        col=adjustcolor("black",alpha.f = 0.1),border = NA)

lines(APP_RED_OUT,type="l",lwd=2,col="red")
polygon(x=c(APP_RED_OUT$YEAR,rev(APP_RED_OUT$YEAR)),y=c(APP_RED_OUT$OBS-APP_RED_OUTsd$OBS,rev(APP_RED_OUT$OBS+APP_RED_OUTsd$OBS)),
        col=adjustcolor("red",alpha.f = 0.1),border = NA)

lines(APP_BLUE_OUT,type="l",lwd=2,col="cadetblue")
polygon(x=c(APP_BLUE_OUT$YEAR,rev(APP_BLUE_OUT$YEAR)),y=c(APP_BLUE_OUT$OBS-APP_BLUE_OUTsd$OBS,rev(APP_BLUE_OUT$OBS+APP_BLUE_OUTsd$OBS)),
        col=adjustcolor("cadetblue",alpha.f = 0.1),border = NA)
dev.off()

# h. Panel E ----
DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

# Load pheno
PHENO = read.table("BIS/DATA/PHENO/PHENO_CLUST4.csv",header=T,sep=",")
PARAM = read.table("BIS/DATA/PHENO/OPTIM_CLUST4.csv",header=T,sep=",")

pdf("BIS5/FIGURES/FIGURE3/PANELS/Figure3_abovepanel.pdf",width = 7,height = 6)
Y=rep(DL.f(PAR = PARAM$cluster3,DOY = DOY.seq),5)
X=c(1:(365*5))
plot(y=Y,x=X,type="l",col="black",xlab="Time",ylab="NDVI",ylim=c(0,1),xlim=c(1,1825),lwd=2)
abline(v=seq(1,10000,365),col="grey",lty=2)
Y=rep(DL.f(PAR = PARAM$cluster2,DOY = DOY.seq),5)
X=c(1:(365*5))
lines(y=Y,x=X,type="l",col="red",lwd=2)

Y=rep(DL.f(PAR = PARAM$cluster4,DOY = DOY.seq),5)
rep()
X=c(1:(365*5))
lines(y=Y,x=X,type="l",col="cadetblue",lwd=2)
dev.off()
# ---- III. FIGURE 4 ----
# a. Prepare data ----
# Load pheno
PHENO = read.table("BIS/DATA/PHENO/PHENO_CLUST4.csv",header=T,sep=",")
PARAM = read.table("BIS/DATA/PHENO/OPTIM_CLUST4.csv",header=T,sep=",")

DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

NDVIyts_CL2 = DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)
NDVIyts_CL3 = DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)
NDVIyts_CL4 = DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)

DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

# Sample NDVI under varying revisit intervals and estimate NDVImax using a retrieving function ---
NSAMP       = 10000         # number of sampling
REVISIT.freq = 4   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 274)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]

# NOBSmin is the minimum number of dates for a given revisit interval  
NOBSmin  <- min(unlist(lapply(DO.tmp,function(x) length(x)))) 

# initialize the matrix to store NDVImax estimates  
NDVImax1.tmp_CL2   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL3   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL4   <- matrix(ncol=NOBSmin,nrow=NSAMP)

DOYgs.tmp_CL2      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL3      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL4      <- matrix(ncol=NOBSmin,nrow=NSAMP)


for (i in 1:NOBSmin){
  print(i)
  
  tmp_CL2     <- sapply(1:NSAMP, function(x) NDVIyts_CL2[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL2 = matrix(ncol=length(tmp_CL2),data=tmp_CL2)}
  tmp.list_CL2 <- split(tmp_CL2, seq(ncol(tmp_CL2)))
  
  tmp_CL3     <- sapply(1:NSAMP, function(x) NDVIyts_CL3[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL3 = matrix(ncol=length(tmp_CL3),data=tmp_CL3)}
  tmp.list_CL3 <- split(tmp_CL3, seq(ncol(tmp_CL3)))
  
  tmp_CL4     <- sapply(1:NSAMP, function(x) NDVIyts_CL4[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL4 = matrix(ncol=length(tmp_CL4),data=tmp_CL4)}
  tmp.list_CL4 <- split(tmp_CL4, seq(ncol(tmp_CL4)))
  
  
  DOYgs.f            <- function(x) length(x)
  DOYgs.tmp_CL2[,i]    <- as.numeric(unlist(lapply(tmp.list_CL2,DOYgs.f)))
  DOYgs.tmp_CL3[,i]    <- as.numeric(unlist(lapply(tmp.list_CL3,DOYgs.f)))
  DOYgs.tmp_CL4[,i]    <- as.numeric(unlist(lapply(tmp.list_CL4,DOYgs.f)))
  
  NDVImax.f1          <- function(x) max(x)
  NDVImax1.tmp_CL2[,i]  <- as.numeric(unlist(lapply(tmp.list_CL2,NDVImax.f1)))
  NDVImax1.tmp_CL3[,i]  <- as.numeric(unlist(lapply(tmp.list_CL3,NDVImax.f1)))
  NDVImax1.tmp_CL4[,i]  <- as.numeric(unlist(lapply(tmp.list_CL4,NDVImax.f1)))

  
}  

NOBSrange   <- 1:10

MAX_CL2        <- matrix(NA,length(NOBSrange))
MAX_CL3        <- matrix(NA,length(NOBSrange))
MAX_CL4        <- matrix(NA,length(NOBSrange))

SD_CL2        <- matrix(NA,length(NOBSrange))
SD_CL3        <- matrix(NA,length(NOBSrange))
SD_CL4        <- matrix(NA,length(NOBSrange))

for(i in 1:length(NOBSrange)){ # number of observations (-1)
  SEL          <- which(DOYgs.tmp_CL2==NOBSrange[i],arr.ind=T)
  DF           <- cbind(SEL,NDVImax1=NDVImax1.tmp_CL2[SEL])
  
  MAX_CL2[i] <- (mean(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  MAX_CL3[i] <- (mean(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  MAX_CL4[i] <- (mean(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
  SD_CL2[i] <- (sd(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  SD_CL3[i] <- (sd(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  SD_CL4[i] <- (sd(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
}

# b. Make figure ----
pdf("BIS5/FIGURES/FIGURE3/PANELS/Figure3_biais.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))

plot(MAX_CL2,ylim=c(80,100),type="o",col="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(1,8),
     ylab="Relative underestimation of NDVImax (%)",xaxt="n",
     xlab="Useable observations")

axis(1,1:15)
abline(v=1:15,col="grey",lty=2)
lines(MAX_CL2,ylim=c(65,100),type="o",col="red",bg="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:10,10:0),y=c(c(40,MAX_CL2-SD_CL2),rep(100,11)),border=NA,col=adjustcolor("red",alpha.f = 0.1))

lines(MAX_CL3,ylim=c(65,100),type="o",col="black",bg="black",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:10,10:0),y=c(c(60,MAX_CL3-SD_CL3),rep(100,11)),border=NA,col=adjustcolor("black",alpha.f = 0.1))

lines(MAX_CL4,ylim=c(65,100),type="o",col="cadetblue",bg="cadetblue",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:10,10:0),y=c(c(20,MAX_CL4-SD_CL4),rep(100,11)),border=NA,col=adjustcolor("cadetblue",alpha.f = 0.1))

abline(h=100,lwd=2)
dev.off()

# ---- IV. FIGURE 5 ----
# a. Prepare data ----
DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/  (1+exp(-PAR[6]*(DOY-PAR[4])))  )
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}
DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

# Load pheno
PHENO = read.table("BIS/DATA/PHENO/PHENO_CLUST4.csv",header=T,sep=",")
PARAM = read.table("BIS/DATA/PHENO/OPTIM_CLUST4.csv",header=T,sep=",")

# Max cluster 2 red : 0.5535702
# Max cluster 3 black : 0.7422723
# Max cluster 4 blue : 0.2029945

# Black - Prepare data
OBS_BLACK_OUT = read.csv("D:/BIS2/DATA/CSV/Observations_Black_Outside.csv",sep=";")
OBS_BLACK_OVL = read.csv("BIS2/DATA/CSV/Observations_Black_Overlap.csv",sep=";")

NDVIyts_CL3 = DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)

NSAMP       = 1000       # number of sampling
REVISIT.freq = 4   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]
FOR.SAMP

#NDVI_BLACK_OUT = OBS_BLACK_OUT;NDVI_BLACK_OUT[] = NA
#NDVI_BLACK_OVL = OBS_BLACK_OVL;NDVI_BLACK_OVL[] = NA

NDVI_BLACK_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_Black_Outside.csv",header=T,sep=";")
NDVI_BLACK_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_Black_Overlap.csv",header=T,sep=";")


for(y in 84045:nrow(OBS_BLACK_OUT)){
  print(y)
  for(i in 3:ncol(OBS_BLACK_OUT)){
    
    # Outside
    if(OBS_BLACK_OUT[y,i] == 0){NDVI_BLACK_OUT[y,i] = NA}else{
      
      tmp_CL3_OUT     <- sapply(1:NSAMP, function(x) NDVIyts_CL3[sample(FOR.SAMP,OBS_BLACK_OUT[y,i])])
      if(OBS_BLACK_OUT[y,i] == 1){tmp_CL3_OUT = matrix(ncol=length(tmp_CL3_OUT),data=tmp_CL3_OUT)}
      tmp.list_CL3_OUT <- split(tmp_CL3_OUT, seq(ncol(tmp_CL3_OUT)))
      
      NDVImax.f1          <- function(x) max(x)
      NDVI_BLACK_OUT[y,i]  <- mean(as.numeric(unlist(lapply(tmp.list_CL3_OUT,NDVImax.f1))))
      
    }
    
    # Overlap
    if(OBS_BLACK_OVL[y,i] == 0){NDVI_BLACK_OVL[y,i] = NA}else{
      tmp_CL3_OVL     <- sapply(1:NSAMP, function(x) NDVIyts_CL3[sample(FOR.SAMP,OBS_BLACK_OVL[y,i])])
      if(OBS_BLACK_OVL[y,i] == 1){tmp_CL3_OVL = matrix(ncol=length(tmp_CL3_OVL),data=tmp_CL3_OVL)}
      tmp.list_CL3_OVL <- split(tmp_CL3_OVL, seq(ncol(tmp_CL3_OVL)))
      
      NDVImax.f1          <- function(x) max(x)
      NDVI_BLACK_OVL[y,i]  <- mean(as.numeric(unlist(lapply(tmp.list_CL3_OVL,NDVImax.f1))))
    }
  }
}

write.table(NDVI_BLACK_OVL,"BIS2/DATA/CSV/FalseNDVImax_Black_Overlap.csv",quote=F,sep=";",row.names=F)
write.table(NDVI_BLACK_OUT,"BIS2/DATA/CSV/FalseNDVImax_Black_Outside.csv",quote=F,sep=";",row.names=F)

# Black - Compute trends
BLACK_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_Black_Outside.csv",header=T,sep=";")
BLACK_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_Black_Overlap.csv",header=T,sep=";")

BLACK_OUT_SLP = data.frame(ID = BLACK_OUT$ID,DEM = BLACK_OUT$DEM,SLP = NA, PVAL = NA)
BLACK_OVL_SLP = data.frame(ID = BLACK_OVL$ID,DEM = BLACK_OVL$DEM,SLP = NA, PVAL = NA)

for(i in 1:nrow(BLACK_OUT)){
  print(i)
  dataf = na.omit(data.frame(NDVI = as.numeric(BLACK_OUT[i,3:40]), YEAR = 1984:2021))
  MOD = mblm(NDVI~YEAR,dataframe = dataf)
  BLACK_OUT_SLP$SLP[i] = coefficients(MOD)[[2]]
  MK = MannKendall(zoo(dataf$NDVI))
  BLACK_OUT_SLP$PVAL[i] = as.numeric(MK$sl)
}

for(i in 1:nrow(BLACK_OVL)){
  print(i)
  dataf = na.omit(data.frame(NDVI = as.numeric(BLACK_OVL[i,3:40]), YEAR = 1984:2021))
  MOD = mblm(NDVI~YEAR,dataframe = dataf)
  BLACK_OVL_SLP$SLP[i] = coefficients(MOD)[[2]]
  MK = MannKendall(zoo(dataf$NDVI))
  BLACK_OVL_SLP$PVAL[i] = as.numeric(MK$sl)
}

write.table(BLACK_OVL_SLP,"BIS2/DATA/CSV/Slope_Black_Overlap.csv",quote=F,sep=";",row.names=F)
write.table(BLACK_OUT_SLP,"BIS2/DATA/CSV/Slope_Black_Outside.csv",quote=F,sep=";",row.names=F)


# Red - Prepare data
OBS_RED_OUT = read.csv("BIS2/DATA/CSV/Observations_Red_Outside.csv",sep=";")
OBS_RED_OVL = read.csv("BIS2/DATA/CSV/Observations_Red_Overlap.csv",sep=";")

NDVIyts_CL2 = DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)

NSAMP       = 1000       # number of sampling
REVISIT.freq = 4   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]
FOR.SAMP

NDVI_RED_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_RED_Outside.csv",header=T,sep=";")
NDVI_RED_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_RED_Overlap.csv",header=T,sep=";")

for(y in 1:nrow(OBS_RED_OUT)){
  print(y)
  for(i in 3:ncol(OBS_RED_OUT)){
    
    # Outside
    if(OBS_RED_OUT[y,i] == 0){NDVI_RED_OUT[y,i] = NA}else{
      tmp_CL2_OUT     <- sapply(1:NSAMP, function(x) NDVIyts_CL2[sample(FOR.SAMP,OBS_RED_OUT[y,i])])
      if(OBS_RED_OUT[y,i] == 1){tmp_CL2_OUT = matrix(ncol=length(tmp_CL2_OUT),data=tmp_CL2_OUT)}
      tmp.list_CL2_OUT <- split(tmp_CL2_OUT, seq(ncol(tmp_CL2_OUT)))
      
      NDVImax.f1          <- function(x) max(x)
      NDVI_RED_OUT[y,i]  <- mean(as.numeric(unlist(lapply(tmp.list_CL2_OUT,NDVImax.f1))))
    }
    
    # Overlap
    if(OBS_RED_OVL[y,i] == 0){NDVI_RED_OVL[y,i] = NA}else{
      tmp_CL2_OVL     <- sapply(1:NSAMP, function(x) NDVIyts_CL2[sample(FOR.SAMP,OBS_RED_OVL[y,i])])
      if(OBS_RED_OVL[y,i] == 1){tmp_CL2_OVL = matrix(ncol=length(tmp_CL2_OVL),data=tmp_CL2_OVL)}
      tmp.list_CL2_OVL <- split(tmp_CL2_OVL, seq(ncol(tmp_CL2_OVL)))
      
      NDVImax.f1          <- function(x) max(x)
      NDVI_RED_OVL[y,i]  <- mean(as.numeric(unlist(lapply(tmp.list_CL2_OVL,NDVImax.f1))))
    }
  }
}

write.table(NDVI_RED_OVL,"BIS2/DATA/CSV/FalseNDVImax_Red_Overlap.csv",quote=F,sep=";",row.names=F)
write.table(NDVI_RED_OUT,"BIS2/DATA/CSV/FalseNDVImax_Red_Outside.csv",quote=F,sep=";",row.names=F)

# Red - Compute trends
RED_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_Red_Outside.csv",header=T,sep=";")
RED_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_Red_Overlap.csv",header=T,sep=";")

RED_OUT_SLP = data.frame(ID = RED_OUT$ID,DEM = RED_OUT$DEM,SLP = NA, PVAL = NA)
RED_OVL_SLP = data.frame(ID = RED_OVL$ID,DEM = RED_OVL$DEM,SLP = NA, PVAL = NA)

for(i in 1:nrow(RED_OUT)){
  print(i)
  dataf = na.omit(data.frame(NDVI = as.numeric(RED_OUT[i,3:40]), YEAR = 1984:2021))
  MOD = mblm(NDVI~YEAR,dataframe = dataf)
  RED_OUT_SLP$SLP[i] = coefficients(MOD)[[2]]
  MK = MannKendall(zoo(dataf$NDVI))
  RED_OUT_SLP$PVAL[i] = as.numeric(MK$sl)
}

for(i in 1:nrow(RED_OVL)){
  print(i)
  dataf = na.omit(data.frame(NDVI = as.numeric(RED_OVL[i,3:40]), YEAR = 1984:2021))
  MOD = mblm(NDVI~YEAR,dataframe = dataf)
  RED_OVL_SLP$SLP[i] = coefficients(MOD)[[2]]
  MK = MannKendall(zoo(dataf$NDVI))
  RED_OVL_SLP$PVAL[i] = as.numeric(MK$sl)
}

write.table(RED_OVL_SLP,"BIS2/DATA/CSV/Slope_Red_Overlap.csv",quote=F,sep=";",row.names=F)
write.table(RED_OUT_SLP,"BIS2/DATA/CSV/Slope_Red_Outside.csv",quote=F,sep=";",row.names=F)


# Blue - Prepare data 
OBS_BLUE_OUT = read.csv("BIS2/DATA/CSV/Observations_BLUE_Outside.csv",sep=";")
OBS_BLUE_OVL = read.csv("BIS2/DATA/CSV/Observations_BLUE_Overlap.csv",sep=";")

NDVIyts_CL4 = DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)

NSAMP       = 1000       # number of sampling
REVISIT.freq = 4   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]
FOR.SAMP

NDVI_BLUE_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_BLUE_Outside.csv",header=T,sep=";")
NDVI_BLUE_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_BLUE_Overlap.csv",header=T,sep=";")

for(y in 35162:nrow(OBS_BLUE_OUT)){
  print(y)
  for(i in 3:ncol(OBS_BLUE_OUT)){
    
    # Outside
    if(OBS_BLUE_OUT[y,i] == 0){NDVI_BLUE_OUT[y,i] = NA}else{
      tmp_CL4_OUT     <- sapply(1:NSAMP, function(x) NDVIyts_CL4[sample(FOR.SAMP,OBS_BLUE_OUT[y,i])])
      if(OBS_BLUE_OUT[y,i] == 1){tmp_CL4_OUT = matrix(ncol=length(tmp_CL4_OUT),data=tmp_CL4_OUT)}
      tmp.list_CL4_OUT <- split(tmp_CL4_OUT, seq(ncol(tmp_CL4_OUT)))
      
      NDVImax.f1          <- function(x) max(x)
      NDVI_BLUE_OUT[y,i]  <- mean(as.numeric(unlist(lapply(tmp.list_CL4_OUT,NDVImax.f1))))
    }
    
    # Overlap
    if(OBS_BLUE_OVL[y,i] == 0){NDVI_BLUE_OVL[y,i] = NA}else{
      tmp_CL4_OVL     <- sapply(1:NSAMP, function(x) NDVIyts_CL4[sample(FOR.SAMP,OBS_BLUE_OVL[y,i])])
      if(OBS_BLUE_OVL[y,i] == 1){tmp_CL4_OVL = matrix(ncol=length(tmp_CL4_OVL),data=tmp_CL4_OVL)}
      tmp.list_CL4_OVL <- split(tmp_CL4_OVL, seq(ncol(tmp_CL4_OVL)))
      
      NDVImax.f1          <- function(x) max(x)
      NDVI_BLUE_OVL[y,i]  <- mean(as.numeric(unlist(lapply(tmp.list_CL4_OVL,NDVImax.f1))))
    }
  }
}

write.table(NDVI_BLUE_OVL,"BIS2/DATA/CSV/FalseNDVImax_BLUE_Overlap.csv",quote=F,sep=";",row.names=F)
write.table(NDVI_BLUE_OUT,"BIS2/DATA/CSV/FalseNDVImax_BLUE_Outside.csv",quote=F,sep=";",row.names=F)

# Blue - Compute trends 
BLUE_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_Blue_Outside.csv",header=T,sep=";")
BLUE_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_Blue_Overlap.csv",header=T,sep=";")

BLUE_OUT_SLP = data.frame(ID = BLUE_OUT$ID,DEM = BLUE_OUT$DEM,SLP = NA, PVAL = NA)
BLUE_OVL_SLP = data.frame(ID = BLUE_OVL$ID,DEM = BLUE_OVL$DEM,SLP = NA, PVAL = NA)

for(i in 1:nrow(BLUE_OUT)){
  print(i)
  dataf = na.omit(data.frame(NDVI = as.numeric(BLUE_OUT[i,3:40]), YEAR = 1984:2021))
  MOD = mblm(NDVI~YEAR,dataframe = dataf)
  BLUE_OUT_SLP$SLP[i] = coefficients(MOD)[[2]]
  MK = MannKendall(zoo(dataf$NDVI))
  BLUE_OUT_SLP$PVAL[i] = as.numeric(MK$sl)
}

for(i in 1:nrow(BLUE_OVL)){
  print(i)
  dataf = na.omit(data.frame(NDVI = as.numeric(BLUE_OVL[i,3:40]), YEAR = 1984:2021))
  MOD = mblm(NDVI~YEAR,dataframe = dataf)
  BLUE_OVL_SLP$SLP[i] = coefficients(MOD)[[2]]
  MK = MannKendall(zoo(dataf$NDVI))
  BLUE_OVL_SLP$PVAL[i] = as.numeric(MK$sl)
}

write.table(BLUE_OVL_SLP,"BIS2/DATA/CSV/Slope_Blue_Overlap.csv",quote=F,sep=";",row.names=F)
write.table(BLUE_OUT_SLP,"BIS2/DATA/CSV/Slope_Blue_Outside.csv",quote=F,sep=";",row.names=F)


# b. Black - Panel 1 ----
OBS_BLACK_OUT = read.csv("BIS2/DATA/CSV/Observations_Black_Outside.csv",sep=";")
OBS_BLACK_OVL = read.csv("BIS2/DATA/CSV/Observations_Black_Overlap.csv",sep=";")
BLACK_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_Black_Outside.csv",header=T,sep=";")
BLACK_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_Black_Overlap.csv",header=T,sep=";")

z=100000

DF_BLACK_OUT = data.frame(MEAN = as.numeric(apply(BLACK_OUT[1:z,3:40],2,function(x){mean(x,na.rm=T)})),
                          SD = as.numeric(apply(BLACK_OUT[1:z,3:40],2,function(x){sd(x,na.rm=T)})),
                          YEAR = 1984:2021)

DF_BLACK_OVL = data.frame(MEAN = as.numeric(apply(BLACK_OVL[1:z,3:40],2,function(x){mean(x,na.rm=T)})),
                          SD = as.numeric(apply(BLACK_OVL[1:z,3:40],2,function(x){sd(x,na.rm=T)})),
                          YEAR = 1984:2021)

pdf("BIS5/FIGURES/FIGURE5/PANELS/Figure5_timeseries_black.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))
plot(y=DF_BLACK_OUT$MEAN,x=DF_BLACK_OUT$YEAR,type="o",ylim=c(0.7,0.744),col="#7b3294",pch=17,cex=0.8,
     ylab="Estimated NDVI max",xlab="Year",xaxs="i")
polygon(x=c(DF_BLACK_OUT$YEAR,rev(DF_BLACK_OUT$YEAR)),y=c(DF_BLACK_OUT$MEAN-DF_BLACK_OUT$SD,rev(DF_BLACK_OUT$MEAN+DF_BLACK_OUT$SD)),
        col=adjustcolor("#7b3294",alpha.f = 0.1),border = NA)

#MODOUT = mblm(MEAN~YEAR,dataframe = DF_BLACK_OUT)
points(y=DF_BLACK_OVL$MEAN,x=DF_BLACK_OVL$YEAR,type="o",ylim=c(0,1),col="#008837",pch=17,cex=0.8)
polygon(x=c(DF_BLACK_OVL$YEAR,rev(DF_BLACK_OVL$YEAR)),y=c(DF_BLACK_OVL$MEAN-DF_BLACK_OVL$SD,rev(DF_BLACK_OVL$MEAN+DF_BLACK_OVL$SD)),
        col=adjustcolor("#008837",alpha.f = 0.1),border = NA)

abline(h=0.7422,col="black",lty=1,lwd=1)
points(y=rep(0.7422,38),x=1984:2021,col="black",pch=17,cex=0.8)

dev.off()

# c. Red - Panel 2 ----
OBS_RED_OUT = read.csv("BIS2/DATA/CSV/Observations_Red_Outside.csv",sep=";")
OBS_RED_OVL = read.csv("BIS2/DATA/CSV/Observations_Red_Overlap.csv",sep=";")
FAKE_RED_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_RED_Outside.csv",header=T,sep=";")
FAKE_RED_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_RED_Overlap.csv",header=T,sep=";")

# Select proper cluster
# RED = 3
# blue = 4
# red = 2
RED_OUT = FAKE_RED_OUT
RED_OVL = FAKE_RED_OVL
z=100000

DF_RED_OUT = data.frame(MEAN = as.numeric(apply(RED_OUT[1:z,3:40],2,function(x){mean(x,na.rm=T)})),
                        SD = as.numeric(apply(RED_OUT[1:z,3:40],2,function(x){sd(x,na.rm=T)})),
                        YEAR = 1984:2021)

DF_RED_OVL = data.frame(MEAN = as.numeric(apply(RED_OVL[1:z,3:40],2,function(x){mean(x,na.rm=T)})),
                        SD = as.numeric(apply(RED_OVL[1:z,3:40],2,function(x){sd(x,na.rm=T)})),
                        YEAR = 1984:2021)


pdf("BIS5/FIGURES/FIGURE5/PANELS/Figure5_timeseries_red.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))
plot(y=DF_RED_OUT$MEAN,x=DF_RED_OUT$YEAR,type="o",ylim=c(0.461,0.558),col="#7b3294",pch=17,cex=0.8,
     ylab="Estimated NDVI max",xlab="Year",xaxs="i")
polygon(x=c(DF_RED_OUT$YEAR,rev(DF_RED_OUT$YEAR)),y=c(DF_RED_OUT$MEAN-DF_RED_OUT$SD,rev(DF_RED_OUT$MEAN+DF_RED_OUT$SD)),
        col=adjustcolor("#7b3294",alpha.f = 0.1),border = NA)
points(y=DF_RED_OVL$MEAN,x=DF_RED_OVL$YEAR,type="o",ylim=c(0,1),col="#008837",pch=17,cex=0.8)
polygon(x=c(DF_RED_OVL$YEAR,rev(DF_RED_OVL$YEAR)),y=c(DF_RED_OVL$MEAN-DF_RED_OVL$SD,rev(DF_RED_OVL$MEAN+DF_RED_OVL$SD)),
        col=adjustcolor("#008837",alpha.f = 0.1),border = NA)

abline(h=0.5535702,col="black",lty=1,lwd=1)
points(y=rep(0.5535702,38),x=1984:2021,col="black",pch=17,cex=0.8)

dev.off()

# d. Blue - Panel 1 ----
OBS_BLUE_OUT = read.csv("BIS2/DATA/CSV/Observations_BLUE_Outside.csv",sep=";")
OBS_BLUE_OVL = read.csv("BIS2/DATA/CSV/Observations_BLUE_Overlap.csv",sep=";")
FAKE_BLUE_OUT = read.table("BIS2/DATA/CSV/FalseNDVImax_BLUE_Outside.csv",header=T,sep=";")
FAKE_BLUE_OVL = read.table("BIS2/DATA/CSV/FalseNDVImax_BLUE_Overlap.csv",header=T,sep=";")

# Select proper cluster
# BLUE = 3
# blue = 4
# BLUE = 2
BLUE_OUT = FAKE_BLUE_OUT
BLUE_OVL = FAKE_BLUE_OVL
z=100000

DF_BLUE_OUT = data.frame(MEAN = as.numeric(apply(BLUE_OUT[1:z,3:40],2,function(x){mean(x,na.rm=T)})),
                         SD = as.numeric(apply(BLUE_OUT[1:z,3:40],2,function(x){sd(x,na.rm=T)})),
                         YEAR = 1984:2021)
DF_BLUE_OVL = data.frame(MEAN = as.numeric(apply(BLUE_OVL[1:z,3:40],2,function(x){mean(x,na.rm=T)})),
                         SD = as.numeric(apply(BLUE_OVL[1:z,3:40],2,function(x){sd(x,na.rm=T)})),
                         YEAR = 1984:2021)


pdf("BIS5/FIGURES/FIGURE5/PANELS/Figure5_timeseries_snowbed.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))
plot(y=DF_BLUE_OUT$MEAN,x=DF_BLUE_OUT$YEAR,type="o",ylim=c(0.125,0.206),col="#7b3294",pch=17,cex=0.8,
     ylab="Estimated NDVI max",xlab="Year",xaxs="i")
polygon(x=c(DF_BLUE_OUT$YEAR,rev(DF_BLUE_OUT$YEAR)),y=c(DF_BLUE_OUT$MEAN-DF_BLUE_OUT$SD,rev(DF_BLUE_OUT$MEAN+DF_BLUE_OUT$SD)),
        col=adjustcolor("#7b3294",alpha.f = 0.1),border = NA)
points(y=DF_BLUE_OVL$MEAN,x=DF_BLUE_OVL$YEAR,type="o",ylim=c(0,1),col="#008837",pch=17,cex=0.8)
polygon(x=c(DF_BLUE_OVL$YEAR,rev(DF_BLUE_OVL$YEAR)),y=c(DF_BLUE_OVL$MEAN-DF_BLUE_OVL$SD,rev(DF_BLUE_OVL$MEAN+DF_BLUE_OVL$SD)),
        col=adjustcolor("#008837",alpha.f = 0.1),border = NA)


abline(h=0.2029945,col="black",lty=1,lwd=1)
points(y=rep(0.2029945,38),x=1984:2021,col="black",pch=17,cex=0.8)
dev.off()

# e. Panel 4,5, 6 ----

BLACK_OUT_SLP = read.table("BIS2/DATA/CSV/Slope_Black_Outside.csv",header=T,sep=";")
BLACK_OVL_SLP = read.table("BIS2/DATA/CSV/Slope_Black_Overlap.csv",header=T,sep=";")

RED_OUT_SLP = read.table("BIS2/DATA/CSV/Slope_Red_Outside.csv",header=T,sep=";")
RED_OVL_SLP = read.table("BIS2/DATA/CSV/Slope_Red_Overlap.csv",header=T,sep=";")

BLUE_OUT_SLP = read.table("BIS2/DATA/CSV/Slope_Blue_Outside.csv",header=T,sep=";")
BLUE_OVL_SLP = read.table("BIS2/DATA/CSV/Slope_Blue_Overlap.csv",header=T,sep=";")

BLUE_OVL_SLPrel = c(c(BLUE_OVL_SLP$SLP/0.2029945)*100)*38
BLUE_OUT_SLPrel = c(c(BLUE_OUT_SLP$SLP/0.2029945)*100)*38

BLACK_OVL_SLPrel = c(c(BLACK_OVL_SLP$SLP/0.7422)*100)*38
BLACK_OUT_SLPrel = c(c(BLACK_OUT_SLP$SLP/0.7422)*100)*38

RED_OVL_SLPrel = c(c(RED_OVL_SLP$SLP/0.5535702)*100)*38
RED_OUT_SLPrel = c(c(RED_OUT_SLP$SLP/0.5535702)*100)*38


pdf("BIS5/FIGURES/FIGURE5/PANELS/Figure5_boxplot_black.pdf",width = 7,height = 6)
par(mfrow=c(1,2))
par(mar=c(5,5,5,0.3))
boxplot(BLACK_OVL_SLP$SLP,BLACK_OUT_SLP$SLP,outline=F,ylim=c(-0.0005,0.002),
        col=c(adjustcolor("#008837",alpha.f = 0.4),adjustcolor("#7b3294",alpha.f = 0.4)),
        border=c(adjustcolor("black",alpha.f = 0.8),adjustcolor("black",alpha.f = 0.8)),xaxt="n")
abline(h=seq(-0.0005,0.0020,0.0005),col="grey70",lty=2)
abline(h=0,lwd=2)
par(mar=c(5,0.3,5,5))

boxplot(BLACK_OVL_SLPrel,BLACK_OUT_SLPrel,outline=F,ylim=c(-7.5,30),
        col=c(adjustcolor("#008837",alpha.f = 0.4),adjustcolor("#7b3294",alpha.f = 0.4)),
        border=c(adjustcolor("black",alpha.f = 0.8),adjustcolor("black",alpha.f = 0.8)),yaxt="n",xaxt="n")
axis(4,at=c(-7.5,0,7.5,15,22.5,30))
abline(h=0,lwd=2)
abline(h=c(-7.5,0,7.5,15,22.5,30),col="grey70",lty=2)
dev.off()


pdf("BIS5/FIGURES/FIGURE5/PANELS/Figure5_boxplot_red.pdf",width = 7,height = 6)
par(mfrow=c(1,2))
par(mar=c(5,5,5,0.3))
boxplot(RED_OVL_SLP$SLP,RED_OUT_SLP$SLP,outline=F,ylim=c(-0.0005,0.002),
        col=c(adjustcolor("#008837",alpha.f = 0.4),adjustcolor("#7b3294",alpha.f = 0.4)),
        border=c(adjustcolor("black",alpha.f = 0.8),adjustcolor("black",alpha.f = 0.8)),xaxt="n")
abline(h=seq(-0.0005,0.0020,0.0005),col="grey70",lty=2)
abline(h=0,lwd=2)
par(mar=c(5,0.3,5,5))

boxplot(RED_OVL_SLPrel,RED_OUT_SLPrel,outline=F,ylim=c(-7.5,30),
        col=c(adjustcolor("#008837",alpha.f = 0.4),adjustcolor("#7b3294",alpha.f = 0.4)),
        border=c(adjustcolor("black",alpha.f = 0.8),adjustcolor("black",alpha.f = 0.8)),yaxt="n",xaxt="n")
axis(4,at=c(-7.5,0,7.5,15,22.5,30))
abline(h=0,lwd=2)
abline(h=c(-7.5,0,7.5,15,22.5,30),col="grey70",lty=2)
dev.off()

pdf("BIS5/FIGURES/FIGURE5/PANELS/Figure5_boxplot_orange.pdf",width = 7,height = 6)
par(mfrow=c(1,2))
par(mar=c(5,5,5,0.3))
boxplot(BLUE_OVL_SLP$SLP,BLUE_OUT_SLP$SLP,outline=F,ylim=c(-0.0005,0.002),
        col=c(adjustcolor("#008837",alpha.f = 0.4),adjustcolor("#7b3294",alpha.f = 0.4)),
        border=c(adjustcolor("black",alpha.f = 0.8),adjustcolor("black",alpha.f = 0.8)),xaxt="n")
abline(h=seq(-0.0005,0.0020,0.0005),col="grey70",lty=2)
abline(h=0,lwd=2)
par(mar=c(5,0.3,5,5))

boxplot(BLUE_OVL_SLPrel,BLUE_OUT_SLPrel,outline=F,ylim=c(-7.5,30),
        col=c(adjustcolor("#008837",alpha.f = 0.4),adjustcolor("#7b3294",alpha.f = 0.4)),
        border=c(adjustcolor("black",alpha.f = 0.8),adjustcolor("black",alpha.f = 0.8)),yaxt="n",xaxt="n")
axis(4,at=c(-7.5,0,7.5,15,22.5,30))
abline(h=0,lwd=2)
abline(h=c(-7.5,0,7.5,15,22.5,30),col="grey70",lty=2)
dev.off()
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# ----------------------------------------------------------------------- #
# ---- V. FIGURE 6 ----
# a. Prepare data ----
# LOAD DATA AND PREPARE TABLE TO COMPUTE ON
OBS_BLACK_OVL = read.csv("BIS2/DATA/CSV/Observations_BLACK_Overlap.csv",sep=";");
OBS_BLACK_OVL=OBS_BLACK_OVL[sample(1:100000,33333),]
OBS_RED_OVL = read.csv("BIS2/DATA/CSV/Observations_RED_Overlap.csv",sep=";");
OBS_RED_OVL=OBS_RED_OVL[sample(1:100000,33333),]
OBS_BLUE_OVL = read.csv("BIS2/DATA/CSV/Observations_BLUE_Overlap.csv",sep=";");
OBS_BLUE_OVL=OBS_BLUE_OVL[sample(1:100000,33334),]

OBSstrat = cbind(rbind(OBS_BLACK_OVL,OBS_RED_OVL,OBS_BLUE_OVL))
OBSstrat$SUM_OBS = apply(OBSstrat[,3:40],1,sum)
hist(OBSstrat$SUM_OBS)
nrow(OBSstrat)

DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}
DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

# CREATE SERIES OF ONSET AND OFFSET

ONSETtmp = rep(seq(131,170,1),2)
ONSET = c(130,ONSETtmp[order(ONSETtmp)])

OFFSETtmp = rep(seq(300,261,-1),2)
OFFSET = c(OFFSETtmp[order(OFFSETtmp,decreasing = T)],260)
GSL = c()

PARAM = c(-0.0199,0.516,150,300,0.076,0.0543)
NDVIyts = DL.f(PARAM,DOY.seq,RESCALE=T)

# b. Figure Sup 2 ----
pdf("BIS5/FIGURES/SUPP/FIG_S2/FIGURE_S2.pdf",width=7,height=6)
plot(NDVIyts,type="n",xlim=c(60,360),xlab="Day of year",ylab="NDVI")


for(i in 1:length(ONSET)){
  ON = ONSET[i]
  OFF = OFFSET[i]
  PARAM = c(-0.0199,0.516,ON,OFF,0.076,0.0543)
  NDVIyts = DL.f(PARAM,DOY.seq,RESCALE=T)
  lines(NDVIyts)
  GSL[i] = OFF - ON
  abline(v=c(152,243))
  
}
dev.off()
print(GSL)

# c. Launch model computation ----
MAT=matrix(ncol = 83,nrow=nrow(OBSstrat))
colnames(MAT) = c("DEM","SUM_OBS",paste0("GSL",seq(170,90,-1)))
MAT = as.data.frame(MAT)
MAT$DEM = OBSstrat$DEM
MAT$SUM_OBS = OBSstrat$SUM_OBS

for(v in 1:length(ONSET)){
  
  # LOAD PARAMETERS
  ON = ONSET[v]
  OFF = OFFSET[v]
  GSL = OFF - ON
  
  # COMPUTE PHENOLOGY
  PARAM = c(-0.0199,0.516,ON,OFF,0.076,0.0543)
  NDVIyts = DL.f(PARAM,DOY.seq,RESCALE=T)
  plot(NDVIyts)
  
  # SAMPLING PARAMETERS
  NSAMP       = 10       # number of sampling
  REVISIT.freq = 4   # revisit interval in days
  
  # SAMPLE DATES
  DO.tmp   <- DO.d(REVISIT.freq)
  for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
  FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]
  
  OBS = as.data.frame(OBSstrat[,3:40])
  MAX = OBS;MAX[] = NA
  
  for(y in 1:nrow(OBS)){
    print(y)
    
    for(z in 1:ncol(OBS)){
      
      if(OBS[y,z] == 0){MAX[y,z] = NA}else{
        
        tmp_CL3_OUT     <- sapply(1:NSAMP, function(x) NDVIyts[sample(FOR.SAMP,OBS[y,z])])
        if(OBS[y,z] == 1){tmp_CL3_OUT = matrix(ncol=length(tmp_CL3_OUT),data=tmp_CL3_OUT)}
        tmp.list_CL3_OUT <- split(tmp_CL3_OUT, seq(ncol(tmp_CL3_OUT)))
        
        NDVImax.f1          <- function(x) max(x)
        MAX[y,z]  <- mean(as.numeric(unlist(lapply(tmp.list_CL3_OUT,NDVImax.f1))))
        
      }
      
      
    }
    dataf = na.omit(data.frame(NDVI = as.numeric(MAX[y,]), YEAR = 1984:2021))
    MAT[y,v+2] = coefficients(mblm(NDVI~YEAR,dataframe = dataf))[[2]]
    
  }
  
  write.table(MAT,paste0("BIS3/DATA/CSV/SimuSensi_GSL",GSL,".csv"),quote=F,sep=";",row.names=F)
  
} # PHENO

# WHAT IS THE TRUE DISTRIBUTION OF GSL AND SUM OBS ?
PIX = fread("BIS2/DATA/PIX/PIX.csv",header=T)
PIXshp = vect(PIX,geom = c("X_ETRS89","Y_ETRS89"), crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs
")
COLD = vect("BIS3/DATA/COLD_REGIONS/CLASSIF/EarthColdRegions.shp")
PIXshp = project(PIXshp,COLD)
PIXshp = mask(PIXshp,COLD)
PIXshp$GSL = PIXshp$OFFSET - PIXshp$ONSET
SUMOBS = rast("BIS2/DATA/OBS/SumOBS.tif")
PIXshp = project(PIXshp,crs(SUMOBS))
SOBS = terra::extract(SUMOBS,PIXshp)[,2]
DF_TRUE = data.frame(GSL = PIXshp$GSL, SOBS = SOBS)
write.table(DF_TRUE,paste0("BIS3/DATA/CSV/Distrib_GSL_OBS.csv"))

# MAKE HISTOGRAM OF TRUE DISTRIBUTION
DF_TRUE = read.table("BIS3/DATA/CSV/Distrib_GSL_OBS.csv")
hist(DF_TRUE$GSL);abline(v=c(170,90))
hist(DF_TRUE$SOBS)


# c. Smoothed heatmap ----
# Load results
RESULTAT = read.table(paste0("BIS3/DATA/CSV/SimuSensi_GSL100.csv"),sep=";",header=T)
DF_TRUE = read.table("BIS2/DATA/CSV/Distrib_GSL_OBS.csv")

# Prepare data frame for heatmap
SEQ = seq(20,360,10)
SEQi = c()
for(i in 1:c(length(SEQ)-1)){SEQi[i] = paste(SEQ[i],SEQ[i+1],sep="-")}
HMdf = data.frame(SUM_OBS = SEQi,ONSET170 = NA)
RES.ONSET = RESULTAT[,3:ncol(RESULTAT)]

for(y in 1:ncol(RES.ONSET)){
  for(i in 1:nrow(HMdf)){
    
    LimInf = SEQ[i]
    LimSup = SEQ[i+1]
    
    HMdf[i,y+1] = median(RES.ONSET[,y][which(RESULTAT$SUM_OBS <= LimSup & RESULTAT$SUM_OBS >= LimInf)],na.rm=T)
    colnames(HMdf)[y+1] = colnames(RES.ONSET)[y]
    
  }
}
DDD = apply(RESULTAT[,-c(1,2)],2,function(x){median(x,na.rm=T)})
DFGSL = data.frame(GSL = substr(names(DDD),4,6),VAL = as.numeric(DDD))

plot(DFGSL)

HMdf$SUM_OBS = seq(25,355,10)
HMdfON = HMdf
dt2 <- HMdfON %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
dt2$rowname = rep(HMdf$SUM_OBS,length(colnames(HMdfON)))
dt2$colname = as.character(as.numeric(substr(dt2$colname,4,7)))
dt2$rowname = factor(dt2$rowname)
dt2$colname = as.factor(dt2$colname)

# Plot heatmap
coul <- c("darkblue" ,"#135AD9" ,"#107CD7" ,"#089BCE" ,"#16B0B3" ,"#58BC8B" ,"#A3BD6A" ,"#DFB951" ,"#FBCD2D" ,"#F9FB0E","yellow","yellow")



dt2 = na.omit(dt2)
dt2$rowname = as.numeric(as.character(dt2$rowname))
dt2$colname = as.numeric(as.character(dt2$colname))
colnames(dt2) = c("x","y","z")

pdf("BIS4/FIGURES/FIGURE4/HEATMAP.pdf",height=5,width=6)
P = ggplot(dt2, aes(x=x, y=y, z=z)) + 
  geom_tile(data=dt2, aes(fill=z)) +
  geom_contour(colour="white",binwidth=0.00025)+
  scale_x_continuous(limits=c(40,300),expand=c(0,0))+
  scale_y_continuous(limits=c(100,150),expand=c(0,0))+
  scale_fill_gradientn(limits = c(0,0.002),colours = coul)+
  theme(legend.position = "right", legend.key.height = unit(2, "cm"))

P
print(P)
dev.off()

# ---- VI. FIGURE 7 and 8 ----
# a. Launch Random Forest ----
RESULTAT = read.table(paste0("BIS3/DATA/CSV/SimuSensi_GSL100.csv"),sep=";",header=T)[,1:73]

new_df <- data.frame(
  OBS = rep(RESULTAT$SUM_OBS,ncol(RESULTAT)-2),
  GSL = rep(seq(170,100,-1),each = nrow(RESULTAT)),
  BF = unlist(RESULTAT[,3:ncol(RESULTAT)]),
  stringsAsFactors = FALSE
)

AGG = aggregate(new_df$BF,list(new_df$OBS,new_df$GSL),median)
colnames(AGG) = c("OBS","GSL","BF")

data.ini = AGG
data.ini = na.omit(data.ini)

## Dispatch 2/3 between training and evaluation
print('dispatching data 2/3...')
train.indx <- sample(1:nrow(data.ini), size = nrow(data.ini))
data.ini.train <- data.ini[train.indx,]
data.ini.eval <- data.ini[-train.indx,]

## Fitting model
control <- trainControl(method="oob", number=20, search="random")
rf.trained <- train(BF~., data=data.ini.train, method="rf", metric="Rsquared", trControl=control, ntree = 500, importance = T)

# b. Apply model ----
load("BIS3/ModelRF.RData")

PIX = fread("BIS2/DATA/PIX/PIX.csv",header=T)
PIXshp = vect(PIX,geom = c("X_ETRS89","Y_ETRS89"), crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs")
PIXshp$GSL = PIXshp$OFFSET - PIXshp$ONSET
SUMOBS = rast("BIS2/DATA/OBS/SumOBS.tif")
PIXshp = project(PIXshp,crs(SUMOBS))
PIXshp$OBS = terra::extract(SUMOBS,PIXshp)[,2]
TILE = rast("BIS2/DATA/TILE/TILE.tif")
PIXshp$TILE = terra::extract(TILE,PIXshp)[,2]

TOMOD = data.frame(ID = PIXshp$ID250falp, GSL = PIXshp$GSL, OBS = PIXshp$OBS)
TOMOD$predicted <- predict.train(rf.trained, newdata = TOMOD[,-1])
PIXshp$BFpredicted = TOMOD$predicted

SLP_RUMPF = rast("BIS4/DATA/NDVI_sen_stable0_reproj.tif")
SLP_RAW = rast("BIS4/DATA/Greening_AlpesEU_reproj.tif")
PIXshp$slp_mod = terra::extract(SLP_MOD,PIXshp)[,2]
PIXshp$slp_rumpf = terra::extract(SLP_RUMPF,PIXshp)[,2]
PIXshp$slp_raw = terra::extract(SLP_RAW,PIXshp)[,2]
PIXshp$slp_raw = PIXshp$slp_raw/10000

SLP_SC = rast("BIS4/DATA/Greening_AlpesEU_SensorCalib_reproj.tif")
PIXshp$slp_sc = terra::extract(SLP_SC,PIXshp)[,2]
PIXshp$slp_sc = PIXshp$slp_sc/10000

# c.Compute statistics ----

PIXf = PIXshp[,c("ELEVavg","DAH","NDVImax","GSL","OBS","BFpredicted","TILE","slp_rumpf","slp_raw","slp_sc")]
PIXfout = PIXf[which(PIXf$TILE == 1),]
PIXfovl = PIXf[which(PIXf$TILE == 2),]

PIXf$ES = (((PIXf$BFpredicted)/(PIXf$NDVImax))*100)*38
PIXf$RI = ((PIXf$BFpredicted)/(PIXf$slp_raw))*100

AGG = aggregate(PIXf$BFpredicted,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),quantile)
PIX_BF = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$GSL,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),quantile)
PIX_GSL = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$OBS,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),quantile)
PIX_OBS = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$NDVImax,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),quantile)
PIX_NDVI = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$ES,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),quantile)
PIX_ES = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$slp_rumpf,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),function(x){quantile(x,na.rm=T)})
PIX_RUMPF = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$slp_raw,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),function(x){quantile(x,na.rm=T)})
PIX_RAW = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$slp_sc,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),function(x){quantile(x,na.rm=T)})
PIX_SC = data.frame(ELEV = seq(1405,2995,10),AGG$x)

AGG = aggregate(PIXf$RI,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),function(x){quantile(x,na.rm=T)})
PIX_RI = data.frame(ELEV = seq(1405,2995,10),AGG$x)

PIXf$diff = PIXf$slp_sc - PIXf$slp_raw
AGG = aggregate(PIXf$diff,list(cut(PIXf$ELEVavg,seq(1400,3000,10))),function(x){quantile(x,na.rm=T)})
PIX_DIFF = data.frame(ELEV = seq(1405,2995,10),AGG$x)


# e. Figure 7 : Panel A ----

pdf("BIS5/FIGURES/FIGURE7/PANELS/Panel_A.pdf",height=7,width=6)

# OBS
plot(1,ylim=c(1405,2995),xlim=c(0,225),
     pch=".",yaxs="i",col=adjustcolor("cadetblue3",alpha.f = 0.1),type="n",xaxt="n",xlab="",ylab="",yaxt="n")
axis(2,at=seq(1200,3200,200))
axis(3,at=seq(0,250,50))
lines(PIX_OBS$X50.,PIX_OBS$ELEV,col="cadetblue3",lwd=2)
polygon(x=c(PIX_OBS$X25.,rev(PIX_OBS$X75.)),y=c(PIX_OBS$ELEV,rev(PIX_OBS$ELEV)),
        col=adjustcolor("cadetblue3",alpha.f = 0.5),border = "cadetblue3")

# GSL
par(new=T)
plot(1,ylim=c(1405,2995),xlim=c(100,225),
     pch=".",yaxs="i",col=adjustcolor("royalblue",alpha.f = 0.1),type="n",yaxt="n",xaxt="n",xlab="",ylab="")
abline(h=seq(1200,3200,200),lty=2,col="grey")
axis(2,at=seq(1200,3200,200))
axis(1,at=seq(100,225,25))
lines(PIX_GSL$X50.,PIX_GSL$ELEV,col="royalblue",lwd=2)
polygon(x=c(PIX_GSL$X25.,rev(PIX_GSL$X75.)),y=c(PIX_GSL$ELEV,rev(PIX_GSL$ELEV)),
        col=adjustcolor("royalblue",alpha.f = 0.5),border = "royalblue")

dev.off()
# d. Figure 7 : Panel B ----

pdf("BIS5/FIGURES/FIGURE7/PANELS/Panel_B.pdf",height=7,width=6)

plot(1,yaxs="i",ylim=c(1405,2995),xlim=c(0,23),
     pch=".",col=adjustcolor("darksalmon",alpha.f = 0.1),type="n",yaxt="n",xaxt="n",xlab="",ylab="")
abline(h=seq(1200,3200,200),lty=2,col="grey")
axis(1,at=seq(0,20,5))
lines(PIX_ES$X50.,PIX_ES$ELEV,col="darksalmon",lwd=2)
polygon(x=c(PIX_ES$X25.,rev(PIX_ES$X75.)),y=c(PIX_RI$ELEV,rev(PIX_RI$ELEV)),
        col=adjustcolor("darksalmon",alpha.f = 0.5),border = "darksalmon")

# GSL
par(new=T)
plot(1,yaxs="i",ylim=c(1405,2995),xlim=c(0,0.0025),
     pch=".",col=adjustcolor("darkolivegreen3",alpha.f = 0.1),type="n",yaxt="n",xaxt="n",xlab="",ylab="")
abline(h=seq(1200,3200,200),lty=2,col="grey")
axis(4,at=seq(1200,3200,200))
axis(3,at=seq(0,0.002,0.0005))
lines(PIX_BF$X50.,PIX_BF$ELEV,col="darkolivegreen3",lwd=2)
polygon(x=c(PIX_BF$X25.,rev(PIX_BF$X75.)),y=c(PIX_BF$ELEV,rev(PIX_BF$ELEV)),
        col=adjustcolor("darkolivegreen3",alpha.f = 0.5),border = "darkolivegreen3")

dev.off()

# e. Figure 8 ----

pdf("BIS5/FIGURES/FIGURE8/Fig8temp.pdf",height=7,width=6)

# OBS
plot(1,ylim=c(1405,2995),xlim=c(0,0.005),
     pch=".",yaxs="i",col=adjustcolor("darkgreen",alpha.f = 0.1),type="n",xaxt="n",xlab="",ylab="",yaxt="n")
axis(4,at=seq(1200,3200,200))
axis(2,at=seq(1200,3200,200))
axis(1,at=seq(0,0.05,0.001))
lines(PIX_RAW$X50.,PIX_RAW$ELEV,col="darkgreen",lwd=2)
polygon(x=c(PIX_RAW$X25.,rev(PIX_RAW$X75.)),y=c(PIX_RAW$ELEV,rev(PIX_RAW$ELEV)),
        col=adjustcolor("darkgreen",alpha.f = 0.5),border = "darkgreen")

# GSL
par(new=T)
plot(1,ylim=c(1405,2995),xlim=c(0,100),
     pch=".",yaxs="i",col=adjustcolor("red",alpha.f = 0.1),type="n",yaxt="n",xaxt="n",xlab="",ylab="")
abline(h=seq(1200,3200,200),lty=2,col="grey")
axis(3,at=seq(0,100,10))
lines(PIX_RI$X50.,PIX_RI$ELEV,col="red",lwd=2)
polygon(x=c(PIX_RI$X25.,rev(PIX_RI$X75.)),y=c(PIX_RI$ELEV,rev(PIX_RI$ELEV)),
        col=adjustcolor("red",alpha.f = 0.5),border = "red")

dev.off()

# ---- XXX. FIGURE EXACT DOY TIME SERIES ----
setwd("F:/")
# Load packages for data handling etc.
library(sf)
library(dplyr)
library(purrr)
library(data.table)
library(stringr)
library(rgee)
library(LandsatTS)
library(googledrive)
library(terra)

# Intialize the Earth Engine with rgee
ee_Initialize()

# Extract data
files_exported_from_EE = list.files("F:/Publications/2023/HIERARCHY/DATA/SIMON/",full.names=T)
lsat.dt <- do.call("rbind", lapply(files_exported_from_EE, fread))
lsat.dt <- lsat_format_data(lsat.dt)

# Clean data using mask
lsat.dt.clean <- lsat_clean_data(lsat.dt, geom.max = 12, cloud.max = 80, sza.max = 60, filter.cfmask.snow = T, filter.cfmask.water = T, filter.jrc.water = T)



# Load pheno
PHENO = read.table("D:/BIS/DATA/PHENO/PHENO_CLUST4.csv",header=T,sep=",")
PARAM = read.table("D:/BIS/DATA/PHENO/OPTIM_CLUST4.csv",header=T,sep=",")

DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}
DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

NDVIyts_CL2 = DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)
NDVIyts_CL3 = DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)
NDVIyts_CL4 = DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)

DOYts = aggregate(lsat.dt.clean$doy,list(lsat.dt.clean$sample.id,lsat.dt.clean$year),c)
DOYtsu = unique(DOYts$x)

pdf("D:/BIS5/FIGURES/SUPP/FIG_S3/hist_obs_cleanobs.pdf",height=5,width=6)
hist(unlist(DOYtsu),main="",xlab="Observations day of year",xaxs="i",yaxs="i",breaks=17)
abline(v=153)
abline(v=182)
abline(v=213)
box()
lines(NDVIyts_CL4*100000,col="blue",lwd=2)
lines(NDVIyts_CL2*60000,col="black",lwd=2)
lines(NDVIyts_CL3*60000,col="red",lwd=2)

# FIGURE 4 - ACTUAL DOY TIME SERIES ----
STOCK_CL2_DOY = data.frame(LENGTH = rep(NA,length(DOYtsu)), NDVImax = NA)
STOCK_CL3_DOY = data.frame(LENGTH = rep(NA,length(DOYtsu)), NDVImax = NA)
STOCK_CL4_DOY = data.frame(LENGTH = rep(NA,length(DOYtsu)), NDVImax = NA)

for(i in 1:length(DOYtsu)){
  print(i)
  
  DOYc = unique(DOYtsu[[i]])
  STOCK_CL2_DOY$LENGTH[i] = length(DOYc)
  NDVIc = NDVIyts_CL2[DOYc]
  STOCK_CL2_DOY$NDVImax[i] = max(NDVIc)
  
  DOYc = unique(DOYtsu[[i]])
  STOCK_CL3_DOY$LENGTH[i] = length(DOYc)
  NDVIc = NDVIyts_CL3[DOYc]
  STOCK_CL3_DOY$NDVImax[i] = max(NDVIc)
  
  DOYc = unique(DOYtsu[[i]])
  STOCK_CL4_DOY$LENGTH[i] = length(DOYc)
  NDVIc = NDVIyts_CL4[DOYc]
  STOCK_CL4_DOY$NDVImax[i] = max(NDVIc)
  
}

STOCK_CL2_DOY$DIFF = (STOCK_CL2_DOY$NDVImax/max(NDVIyts_CL2))*100
STOCK_CL3_DOY$DIFF = (STOCK_CL3_DOY$NDVImax/max(NDVIyts_CL3))*100
STOCK_CL4_DOY$DIFF = (STOCK_CL4_DOY$NDVImax/max(NDVIyts_CL4))*100

MAX_CL2_DOY = aggregate(STOCK_CL2_DOY$DIFF,list(STOCK_CL2_DOY$LENGTH),mean)
MAX_CL3_DOY = aggregate(STOCK_CL3_DOY$DIFF,list(STOCK_CL3_DOY$LENGTH),mean)
MAX_CL4_DOY = aggregate(STOCK_CL4_DOY$DIFF,list(STOCK_CL4_DOY$LENGTH),mean)


# FIGURE 4 - RANDOMIZATION OF DOY ----
STOCK_CL2_RAND = data.frame(LENGTH = rep(NA,length(DOYtsu)), NDVImax = NA)
STOCK_CL3_RAND = data.frame(LENGTH = rep(NA,length(DOYtsu)), NDVImax = NA)
STOCK_CL4_RAND = data.frame(LENGTH = rep(NA,length(DOYtsu)), NDVImax = NA)

DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

# Sample NDVI under varying revisit intervals and estimate NDVImax using a retrieving function ---
REVISIT.freq = 4   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]

for(i in 1:length(DOYtsu)){
  print(i)
  
  DOYc = unique(DOYtsu[[i]])
  LEN = length(DOYc)
  STOCK_CL2_RAND$LENGTH[i] = LEN
  STOCK_CL2_RAND$NDVImax[i] = max(NDVIyts_CL2[sample(FOR.SAMP,LEN)])
  
  DOYc = unique(DOYtsu[[i]])
  LEN = length(DOYc)
  STOCK_CL3_RAND$LENGTH[i] = LEN
  STOCK_CL3_RAND$NDVImax[i] = max(NDVIyts_CL3[sample(FOR.SAMP,LEN)])
  
  DOYc = unique(DOYtsu[[i]])
  LEN = length(DOYc)
  STOCK_CL4_RAND$LENGTH[i] = LEN
  STOCK_CL4_RAND$NDVImax[i] = max(NDVIyts_CL4[sample(FOR.SAMP,LEN)])
  
}

STOCK_CL2_RAND$DIFF = (STOCK_CL2_RAND$NDVImax/max(NDVIyts_CL2))*100
STOCK_CL3_RAND$DIFF = (STOCK_CL3_RAND$NDVImax/max(NDVIyts_CL3))*100
STOCK_CL4_RAND$DIFF = (STOCK_CL4_RAND$NDVImax/max(NDVIyts_CL4))*100

MAX_CL2_RAND = aggregate(STOCK_CL2_RAND$DIFF,list(STOCK_CL2_RAND$LENGTH),mean)
MAX_CL3_RAND = aggregate(STOCK_CL3_RAND$DIFF,list(STOCK_CL3_RAND$LENGTH),mean)
MAX_CL4_RAND = aggregate(STOCK_CL4_RAND$DIFF,list(STOCK_CL4_RAND$LENGTH),mean)

pdf("BIS5/FIGURES/SUPP/FIG_S3/diff_rand_doy.pdf",width=6,height = 6)
plot(MAX_CL2_DOY,ylim=c(65,100),type="o",col="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(1,8),
     ylab="Relative underestimation of NDVImax (%)",xaxt="n",
     xlab="Useable observations")

axis(1,1:15)
abline(v=1:15,col="grey",lty=2)
lines(MAX_CL2_DOY,ylim=c(65,100),type="o",col="red",bg="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))

lines(MAX_CL3_DOY,ylim=c(65,100),type="o",col="black",bg="black",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))

lines(MAX_CL4_DOY,ylim=c(65,100),type="o",col="cadetblue",bg="cadetblue",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))

lines(MAX_CL2_RAND,ylim=c(65,100),type="o",col="red",bg="red",lwd=2,pch=23,cex=2,lty=1,xlim=c(0.5,10))

lines(MAX_CL3_RAND,ylim=c(65,100),type="o",col="black",bg="black",lwd=2,pch=23,cex=2,lty=1,xlim=c(0.5,10))

lines(MAX_CL4_RAND,ylim=c(65,100),type="o",col="cadetblue",bg="cadetblue",lwd=2,pch=23,cex=2,lty=1,xlim=c(0.5,10))

legend(x=3,y=75,lty=c(1,2),legend=c("With randomization","Actual DOY time series"),cex=1.1,lwd=2)

dev.off()





boxplot(STOCK_CL2_RAND$DIFF~STOCK_CL2_RAND$LENGTH,outline=F,ylim=c(60,100))
boxplot(STOCK_CL2_DOY$DIFF~STOCK_CL2_DOY$LENGTH,outline=F,ylim=c(60,100))

boxplot(STOCK_CL3_RAND$DIFF~STOCK_CL3_RAND$LENGTH,outline=F,ylim=c(60,100))
boxplot(STOCK_CL3_DOY$DIFF~STOCK_CL3_DOY$LENGTH,outline=F,ylim=c(60,100))

par(mfrow = c(1,2))
par(mar=c(2,2,2,2))
boxplot(STOCK_CL4_RAND$DIFF~STOCK_CL4_RAND$LENGTH,outline=F,ylim=c(60,100),xlim=c(0,10))
boxplot(STOCK_CL4_DOY$DIFF~STOCK_CL4_DOY$LENGTH,outline=F,ylim=c(60,100),xlim=c(0,10))

IQR_RAND = as.data.frame(aggregate(STOCK_CL3_RAND$DIFF,list(STOCK_CL3_RAND$LENGTH),quantile)$x); colnames(IQR_RAND) = c("P0","P25","P50","P75","P100")
IQR_DOY = as.data.frame(aggregate(STOCK_CL3_DOY$DIFF,list(STOCK_CL3_DOY$LENGTH),quantile)$x); colnames(IQR_DOY) = c("P0","P25","P50","P75","P100")
plot(IQR_RAND$P50,IQR_DOY$P50,pch=17,xlim=c(95,100),ylim=c(95,100),cex=2,
     xlab="Relative underestimation of NDVImax when DOY randomized",
     ylab="Relative underestimation of NDVImax with true DOY")
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P75, y0 = IQR_DOY$P50, y1 = IQR_DOY$P50, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P25, y0 = IQR_DOY$P50, y1 = IQR_DOY$P50, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P50, y0 = IQR_DOY$P50, y1 = IQR_DOY$P75, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P50, y0 = IQR_DOY$P50, y1 = IQR_DOY$P25, angle=90,length = 0.1,lty=2)
abline(0,1)

IQR_RAND = as.data.frame(aggregate(STOCK_CL2_RAND$DIFF,list(STOCK_CL2_RAND$LENGTH),quantile)$x); colnames(IQR_RAND) = c("P0","P25","P50","P75","P100")
IQR_DOY = as.data.frame(aggregate(STOCK_CL2_DOY$DIFF,list(STOCK_CL2_DOY$LENGTH),quantile)$x); colnames(IQR_DOY) = c("P0","P25","P50","P75","P100")
plot(IQR_RAND$P50,IQR_DOY$P50,pch=17,xlim=c(85,100),ylim=c(85,100),cex=2,
     xlab="Relative underestimation of NDVImax when DOY randomized",
     ylab="Relative underestimation of NDVImax with true DOY")
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P75, y0 = IQR_DOY$P50, y1 = IQR_DOY$P50, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P25, y0 = IQR_DOY$P50, y1 = IQR_DOY$P50, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P50, y0 = IQR_DOY$P50, y1 = IQR_DOY$P75, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P50, y0 = IQR_DOY$P50, y1 = IQR_DOY$P25, angle=90,length = 0.1,lty=2)
abline(0,1)

IQR_RAND = as.data.frame(aggregate(STOCK_CL4_RAND$DIFF,list(STOCK_CL4_RAND$LENGTH),quantile)$x); colnames(IQR_RAND) = c("P0","P25","P50","P75","P100")
IQR_DOY = as.data.frame(aggregate(STOCK_CL4_DOY$DIFF,list(STOCK_CL4_DOY$LENGTH),quantile)$x); colnames(IQR_DOY) = c("P0","P25","P50","P75","P100")
plot(IQR_RAND$P50,IQR_DOY$P50,pch=17,xlim=c(60,100),ylim=c(60,100),cex=2,
     xlab="Relative underestimation of NDVImax when DOY randomized",
     ylab="Relative underestimation of NDVImax with true DOY")
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P75, y0 = IQR_DOY$P50, y1 = IQR_DOY$P50, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P25, y0 = IQR_DOY$P50, y1 = IQR_DOY$P50, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P50, y0 = IQR_DOY$P50, y1 = IQR_DOY$P75, angle=90,length = 0.1,lty=2)
arrows(x0 = IQR_RAND$P50, x1 = IQR_RAND$P50, y0 = IQR_DOY$P50, y1 = IQR_DOY$P25, angle=90,length = 0.1,lty=2)
abline(0,1)

# FIGURE 5 - ACTUAL DOY TIME SERIES ----
NDVImax_thermic_doy = data.frame(sample.id = unique(DOYts$Group.1));NDVImax_thermic_doy[,2:39] = NA;colnames(NDVImax_thermic_doy)[2:39] = as.character(1984:2021)
NDVImax_inter_doy = data.frame(sample.id = unique(DOYts$Group.1));NDVImax_inter_doy[,2:39] = NA;colnames(NDVImax_inter_doy)[2:39] = as.character(1984:2021)
NDVImax_snowbed_doy = data.frame(sample.id = unique(DOYts$Group.1));NDVImax_snowbed_doy[,2:39] = NA;colnames(NDVImax_snowbed_doy)[2:39] = as.character(1984:2021)
OBS_thermic_doy = data.frame(sample.id = unique(DOYts$Group.1));OBS_thermic_doy[,2:39] = NA;colnames(OBS_thermic_doy)[2:39] = as.character(1984:2021)
OBS_inter_doy = data.frame(sample.id = unique(DOYts$Group.1));OBS_inter_doy[,2:39] = NA;colnames(OBS_inter_doy)[2:39] = as.character(1984:2021)
OBS_snowbed_doy = data.frame(sample.id = unique(DOYts$Group.1));OBS_snowbed_doy[,2:39] = NA;colnames(OBS_snowbed_doy)[2:39] = as.character(1984:2021)

for(i in 25600:length(unique(DOYts$Group.1))){
  print(i)
  
  # Load doy time series
  TSsamp = DOYts[which(DOYts$Group.1 == unique(DOYts$Group.1)[i]),];colnames(TSsamp) = c("sample.id","year","doy")
  if(any(TSsamp$year == 2022 | TSsamp$year == 2023)){
    TSsamp = TSsamp[-which(TSsamp$year == 2022 | TSsamp$year == 2023),]
  }
  TSsamp$ndvimax = NA
  
  # Add row if years are missing
  TOADD = c(1984:2021)[c(1984:2021) %!in% TSsamp$year]
  
  if(length(TSsamp$year) != 38){
    TOADDdf = data.frame(sample.id = TSsamp$sample.id[1], year = TOADD, doy = NA, ndvimax = NA)
    TSsamp = rbind(TSsamp,TOADDdf)
    TSsamp = TSsamp[order(TSsamp$year),]
  }
  
  # Prepare DF to store
  TSsamp_thermic_doy = TSsamp
  TSsamp_inter_doy = TSsamp
  TSsamp_snowbed_doy = TSsamp
  
  # Obtain annual ndvi max time series
  for(y in 1:nrow(TSsamp)){
    if(any(!is.na(TSsamp$doy[[y]]))){
      TSsamp_thermic_doy$ndvimax[y] = max(NDVIyts_CL3[TSsamp$doy[[y]]])
      TSsamp_inter_doy$ndvimax[y] = max(NDVIyts_CL2[TSsamp$doy[[y]]])
      TSsamp_snowbed_doy$ndvimax[y] = max(NDVIyts_CL4[TSsamp$doy[[y]]])
      
      LEN = length(TSsamp$doy[[y]])
      TSsamp_thermic_doy$number[y] = LEN
      TSsamp_inter_doy$number[y] = LEN
      TSsamp_snowbed_doy$number[y] = LEN
    } else {
      TSsamp_thermic_doy$ndvimax[y] = NA
      TSsamp_inter_doy$ndvimax[y] = NA
      TSsamp_snowbed_doy$ndvimax[y] = NA
      
      TSsamp_thermic_doy$number[y] = NA
      TSsamp_inter_doy$number[y] = NA
      TSsamp_snowbed_doy$number[y] = NA
    }
  }
  
  # Store in new data.frame
  NDVImax_thermic_doy[i,2:39] = TSsamp_thermic_doy$ndvimax
  NDVImax_inter_doy[i,2:39] = TSsamp_inter_doy$ndvimax
  NDVImax_snowbed_doy[i,2:39] = TSsamp_snowbed_doy$ndvimax
  
  OBS_thermic_doy[i,2:39] = TSsamp_thermic_doy$number
  OBS_inter_doy[i,2:39] = TSsamp_inter_doy$number
  OBS_snowbed_doy[i,2:39] = TSsamp_snowbed_doy$number
  
}

# FIGURE 5 - RANDOMIZING DOY TIME SERIES ----
NDVImax_thermic_rand = data.frame(sample.id = unique(DOYts$Group.1));NDVImax_thermic_rand[,2:39] = NA;colnames(NDVImax_thermic_rand)[2:39] = as.character(1984:2021)
NDVImax_inter_rand = data.frame(sample.id = unique(DOYts$Group.1));NDVImax_inter_rand[,2:39] = NA;colnames(NDVImax_inter_rand)[2:39] = as.character(1984:2021)
NDVImax_snowbed_rand = data.frame(sample.id = unique(DOYts$Group.1));NDVImax_snowbed_rand[,2:39] = NA;colnames(NDVImax_snowbed_rand)[2:39] = as.character(1984:2021)

for(i in 1:length(unique(DOYts$Group.1))){
  print(i)
  
  # Load doy time series
  TSsamp = DOYts[which(DOYts$Group.1 == unique(DOYts$Group.1)[i]),];colnames(TSsamp) = c("sample.id","year","doy")
  if(any(TSsamp$year == 2022 | TSsamp$year == 2023)){
    TSsamp = TSsamp[-which(TSsamp$year == 2022 | TSsamp$year == 2023),]
  }
  TSsamp$ndvimax = NA
  
  # Add row if years are missing
  TOADD = c(1984:2021)[c(1984:2021) %!in% TSsamp$year]
  
  if(length(TSsamp$year) != 38){
    TOADDdf = data.frame(sample.id = TSsamp$sample.id[1], year = TOADD, doy = NA, ndvimax = NA)
    TSsamp = rbind(TSsamp,TOADDdf)
    TSsamp = TSsamp[order(TSsamp$year),]
  }
  
  # Prepare DF to store
  TSsamp_thermic_rand = TSsamp
  TSsamp_inter_rand = TSsamp
  TSsamp_snowbed_rand = TSsamp
  
  # Obtain annual ndvi max time series
  for(y in 1:nrow(TSsamp)){
    LEN = length(TSsamp$doy[[y]])
    if(any(!is.na(TSsamp$doy[[y]])) & LEN < 20){
      
      TSsamp_thermic_rand$ndvimax[y] = max(NDVIyts_CL3[sample(FOR.SAMP,LEN)])
      TSsamp_inter_rand$ndvimax[y] = max(NDVIyts_CL2[sample(FOR.SAMP,LEN)])
      TSsamp_snowbed_rand$ndvimax[y] = max(NDVIyts_CL4[sample(FOR.SAMP,LEN)])
    } else {
      TSsamp_thermic_rand$ndvimax[y] = NA
      TSsamp_inter_rand$ndvimax[y] = NA
      TSsamp_snowbed_rand$ndvimax[y] = NA
    }
  }
  
  # Store in new data.frame
  NDVImax_thermic_rand[i,2:39] = TSsamp_thermic_rand$ndvimax
  NDVImax_inter_rand[i,2:39] = TSsamp_inter_rand$ndvimax
  NDVImax_snowbed_rand[i,2:39] = TSsamp_snowbed_rand$ndvimax
  
}




plot(apply(NDVImax_snowbed_doy,2,function(x){mean(x,na.rm=T)})[-1],type="l",ylim=c(0.17,0.21))

abline(h=max(NDVIyts_CL4))
plot(apply(OBS_snowbed_doy,2,function(x){mean(x,na.rm=T)})[-1],type="l")

boxplot(NDVImax_snowbed_doy$'2012'~OBS_snowbed_doy$'2012',outline=F)

plot(apply(NDVImax_snowbed_rand,2,function(x){mean(x,na.rm=T)})[-1])






# ---- XXX. FIGURE SUP ? BIAS REVISIT CHANGE ----
DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year

NDVIyts_CL2 = DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)
NDVIyts_CL3 = DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)
NDVIyts_CL4 = DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)

DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

# Sample NDVI under varying revisit intervals and estimate NDVImax using a retrieving function ---
NSAMP       = 10000         # number of sampling
REVISIT.freq = 16   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]

# NOBSmin is the minimum number of dates for a given revisit interval  
NOBSmin  <- min(unlist(lapply(DO.tmp,function(x) length(x)))) 

# initialize the matrix to store NDVImax estimates  
NDVImax1.tmp_CL2   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL3   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL4   <- matrix(ncol=NOBSmin,nrow=NSAMP)

DOYgs.tmp_CL2      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL3      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL4      <- matrix(ncol=NOBSmin,nrow=NSAMP)


for (i in 1:NOBSmin){
  print(i)
  
  tmp_CL2     <- sapply(1:NSAMP, function(x) NDVIyts_CL2[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL2 = matrix(ncol=length(tmp_CL2),data=tmp_CL2)}
  tmp.list_CL2 <- split(tmp_CL2, seq(ncol(tmp_CL2)))
  
  tmp_CL3     <- sapply(1:NSAMP, function(x) NDVIyts_CL3[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL3 = matrix(ncol=length(tmp_CL3),data=tmp_CL3)}
  tmp.list_CL3 <- split(tmp_CL3, seq(ncol(tmp_CL3)))
  
  tmp_CL4     <- sapply(1:NSAMP, function(x) NDVIyts_CL4[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL4 = matrix(ncol=length(tmp_CL4),data=tmp_CL4)}
  tmp.list_CL4 <- split(tmp_CL4, seq(ncol(tmp_CL4)))
  
  
  DOYgs.f            <- function(x) length(x)
  DOYgs.tmp_CL2[,i]    <- as.numeric(unlist(lapply(tmp.list_CL2,DOYgs.f)))
  DOYgs.tmp_CL3[,i]    <- as.numeric(unlist(lapply(tmp.list_CL3,DOYgs.f)))
  DOYgs.tmp_CL4[,i]    <- as.numeric(unlist(lapply(tmp.list_CL4,DOYgs.f)))
  
  NDVImax.f1          <- function(x) max(x)
  NDVImax1.tmp_CL2[,i]  <- as.numeric(unlist(lapply(tmp.list_CL2,NDVImax.f1)))
  NDVImax1.tmp_CL3[,i]  <- as.numeric(unlist(lapply(tmp.list_CL3,NDVImax.f1)))
  NDVImax1.tmp_CL4[,i]  <- as.numeric(unlist(lapply(tmp.list_CL4,NDVImax.f1)))
  
  
}  

NOBSrange   <- 1:10

MAX_CL2        <- matrix(NA,length(NOBSrange))
MAX_CL3        <- matrix(NA,length(NOBSrange))
MAX_CL4        <- matrix(NA,length(NOBSrange))

SD_CL2        <- matrix(NA,length(NOBSrange))
SD_CL3        <- matrix(NA,length(NOBSrange))
SD_CL4        <- matrix(NA,length(NOBSrange))

for(i in 1:length(NOBSrange)){ # number of observations (-1)
  SEL          <- which(DOYgs.tmp_CL2==NOBSrange[i],arr.ind=T)
  DF           <- cbind(SEL,NDVImax1=NDVImax1.tmp_CL2[SEL])
  
  MAX_CL2[i] <- (mean(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  MAX_CL3[i] <- (mean(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  MAX_CL4[i] <- (mean(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
  SD_CL2[i] <- (sd(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  SD_CL3[i] <- (sd(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  SD_CL4[i] <- (sd(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
}

pdf("BIS5/FIGURES/FIGURE3/PANELS/Figure3_biais.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))

plot(MAX_CL2,ylim=c(60,100),type="o",col="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(1,8),
     ylab="Relative underestimation of NDVImax (%)",xaxt="n",
     xlab="Useable observations")

axis(1,1:15)
abline(v=1:15,col="grey",lty=2)
lines(MAX_CL2,ylim=c(65,100),type="o",col="red",bg="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:10,10:0),y=c(c(40,MAX_CL2-SD_CL2),rep(100,11)),border=NA,col=adjustcolor("red",alpha.f = 0.1))

lines(MAX_CL3,ylim=c(65,100),type="o",col="black",bg="black",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:10,10:0),y=c(c(60,MAX_CL3-SD_CL3),rep(100,11)),border=NA,col=adjustcolor("black",alpha.f = 0.1))

lines(MAX_CL4,ylim=c(65,100),type="o",col="cadetblue",bg="cadetblue",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:10,10:0),y=c(c(20,MAX_CL4-SD_CL4),rep(100,11)),border=NA,col=adjustcolor("cadetblue",alpha.f = 0.1))

abline(h=100,lwd=2)
dev.off()
# ---- XXX. FIGURE SUP 4 AND 5 QUANTILE BIAS ----
DL.f <- function(PAR, DOY,RESCALE=F){
  NDVI  <- PAR[1] + (PAR[2] - PAR[1])*(1/(1+exp(-PAR[5]*(DOY-PAR[3])))-1/(1+exp(-PAR[6]*(DOY-PAR[4]))))
  if(RESCALE) NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI)) # rescaling
  return(NDVI)
}

DOY.seq     = seq(1,365,1)     # Day of (non bissextile) year
# Load pheno
PHENO = read.table("BIS/DATA/PHENO/PHENO_CLUST4.csv",header=T,sep=",")
PARAM = read.table("BIS/DATA/PHENO/OPTIM_CLUST4.csv",header=T,sep=",")

NDVIyts_CL2 = DL.f(PAR = PARAM$cluster2,DOY = DOY.seq)
NDVIyts_CL3 = DL.f(PAR = PARAM$cluster3,DOY = DOY.seq)
NDVIyts_CL4 = DL.f(PAR = PARAM$cluster4,DOY = DOY.seq)

DO.d <- function(REVISIT){
  SAMP <- list()
  for (i in 1:(REVISIT-1)) SAMP[[i]] <- seq(i,365,REVISIT)
  return(SAMP)  # a list of length REVISIT-1 with vectors of observation dates in DOY
}

# Sample NDVI under varying revisit intervals and estimate NDVImax using a retrieving function ---
NSAMP       = 10000         # number of sampling
REVISIT.freq = 4   # revisit interval in days

# sample dates
DO.tmp   <- DO.d(REVISIT.freq)
for(i in 1:length(DO.tmp)){DO.tmp[[i]] = DO.tmp[[i]][which(DO.tmp[[i]] > 152 & DO.tmp[[i]] < 243)]}
FOR.SAMP = DO.tmp[[sample(1:length(DO.tmp),1)]]

# NOBSmin is the minimum number of dates for a given revisit interval  
NOBSmin  <- min(unlist(lapply(DO.tmp,function(x) length(x)))) 

# initialize the matrix to store NDVImax estimates  
NDVImax1.tmp_CL2   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL3   <- matrix(ncol=NOBSmin,nrow=NSAMP)
NDVImax1.tmp_CL4   <- matrix(ncol=NOBSmin,nrow=NSAMP)

DOYgs.tmp_CL2      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL3      <- matrix(ncol=NOBSmin,nrow=NSAMP)
DOYgs.tmp_CL4      <- matrix(ncol=NOBSmin,nrow=NSAMP)


for (i in 1:NOBSmin){
  print(i)
  
  tmp_CL2     <- sapply(1:NSAMP, function(x) NDVIyts_CL2[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL2 = matrix(ncol=length(tmp_CL2),data=tmp_CL2)}
  tmp.list_CL2 <- split(tmp_CL2, seq(ncol(tmp_CL2)))
  
  tmp_CL3     <- sapply(1:NSAMP, function(x) NDVIyts_CL3[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL3 = matrix(ncol=length(tmp_CL3),data=tmp_CL3)}
  tmp.list_CL3 <- split(tmp_CL3, seq(ncol(tmp_CL3)))
  
  tmp_CL4     <- sapply(1:NSAMP, function(x) NDVIyts_CL4[sample(FOR.SAMP,i)])
  if(i == 1){tmp_CL4 = matrix(ncol=length(tmp_CL4),data=tmp_CL4)}
  tmp.list_CL4 <- split(tmp_CL4, seq(ncol(tmp_CL4)))
  
  
  DOYgs.f            <- function(x) length(x)
  DOYgs.tmp_CL2[,i]    <- as.numeric(unlist(lapply(tmp.list_CL2,DOYgs.f)))
  DOYgs.tmp_CL3[,i]    <- as.numeric(unlist(lapply(tmp.list_CL3,DOYgs.f)))
  DOYgs.tmp_CL4[,i]    <- as.numeric(unlist(lapply(tmp.list_CL4,DOYgs.f)))
  
  NDVImax.f1          <- function(x) quantile(x,probs=0.75,type=9)
  NDVImax1.tmp_CL2[,i]  <- as.numeric(unlist(lapply(tmp.list_CL2,NDVImax.f1)))
  NDVImax1.tmp_CL3[,i]  <- as.numeric(unlist(lapply(tmp.list_CL3,NDVImax.f1)))
  NDVImax1.tmp_CL4[,i]  <- as.numeric(unlist(lapply(tmp.list_CL4,NDVImax.f1)))
  
  
}  

NOBSrange   <- 1:15

MAX_CL2        <- matrix(NA,length(NOBSrange))
MAX_CL3        <- matrix(NA,length(NOBSrange))
MAX_CL4        <- matrix(NA,length(NOBSrange))

SD_CL2        <- matrix(NA,length(NOBSrange))
SD_CL3        <- matrix(NA,length(NOBSrange))
SD_CL4        <- matrix(NA,length(NOBSrange))

for(i in 1:length(NOBSrange)){ # number of observations (-1)
  SEL          <- which(DOYgs.tmp_CL2==NOBSrange[i],arr.ind=T)
  DF           <- cbind(SEL,NDVImax1=NDVImax1.tmp_CL2[SEL])
  
  
  MAX_CL2[i] <- (mean(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  MAX_CL3[i] <- (mean(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  MAX_CL4[i] <- (mean(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
  SD_CL2[i] <- (sd(NDVImax1.tmp_CL2[SEL])/max(NDVIyts_CL2))*100
  SD_CL3[i] <- (sd(NDVImax1.tmp_CL3[SEL])/max(NDVIyts_CL3))*100
  SD_CL4[i] <- (sd(NDVImax1.tmp_CL4[SEL])/max(NDVIyts_CL4))*100
  
}

pdf("BIS5/FIGURES/SUPP/FIG_S4/Fig_S4_q75_type9.pdf",width = 7,height = 6)
par(mar=c(5,5,5,5))

plot(MAX_CL2,ylim=c(60,100),type="o",col="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(1,15),
     ylab="Relative underestimation of NDVImax (%)",xaxt="n",
     xlab="Useable observations",main="Quantile 0.75 (Type = 9)")

axis(1,1:15)
abline(v=1:15,col="grey",lty=2)
lines(MAX_CL2,ylim=c(65,200),type="o",col="red",bg="red",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:15,15:0),y=c(c(40,MAX_CL2-SD_CL2),rep(100,16)),border=NA,col=adjustcolor("red",alpha.f = 0.1))

lines(MAX_CL3,ylim=c(65,100),type="o",col="black",bg="black",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:15,15:0),y=c(c(60,MAX_CL3-SD_CL3),rep(100,16)),border=NA,col=adjustcolor("black",alpha.f = 0.1))

lines(MAX_CL4,ylim=c(65,100),type="o",col="cadetblue",bg="cadetblue",lwd=2,pch=23,cex=2,lty=2,xlim=c(0.5,10))
polygon(x=c(0:15,15:0),y=c(c(20,MAX_CL4-SD_CL4),rep(100,16)),border=NA,col=adjustcolor("cadetblue",alpha.f = 0.1))

abline(h=100,lwd=2)
dev.off()
# ---- XXX. FIGURE TEST CROSS CALIBRATION ----
plot(1,ylim=c(1405,2995),xlim=c(-0.001,0.005),
     pch=".",yaxs="i",col=adjustcolor("darkgreen",alpha.f = 0.1),type="n",xlab="Greening trends",ylab="Elevation")
abline(v=0,lwd=2)
lines(PIX_SC$X50.,PIX_SC$ELEV,col="darkgreen",lwd=2)
polygon(x=c(PIX_SC$X25.,rev(PIX_SC$X75.)),y=c(PIX_SC$ELEV,rev(PIX_SC$ELEV)),
        col=adjustcolor("darkgreen",alpha.f = 0.5),border = "darkgreen")
par(new=T)
plot(1,ylim=c(1405,2995),xlim=c(-0.001,0.005),
     pch=".",yaxs="i",col=adjustcolor("orange",alpha.f = 0.1),type="n",xaxt="n",xlab="",ylab="",yaxt="n")
axis(1,at=seq(0,0.05,0.001))
lines(PIX_RAW$X50.,PIX_RAW$ELEV,col="orange",lwd=2)
polygon(x=c(PIX_RAW$X25.,rev(PIX_RAW$X75.)),y=c(PIX_RAW$ELEV,rev(PIX_RAW$ELEV)),
        col=adjustcolor("orange",alpha.f = 0.5),border = "orange")

plot(PIXf$diff,PIXf$slp_raw,xlim=c(-0.004,0.001),ylim=c(-0.001,0.005),pch=".")
abline(v=-0.001,col="red")
plot(1,ylim=c(1405,2995),xlim=c(-0.002,0.002),
     pch=".",yaxs="i",col=adjustcolor("darkgreen",alpha.f = 0.1),type="n",xaxt="n",xlab="",ylab="",yaxt="n")
axis(4,at=seq(1200,3200,200))
axis(1,at=seq(-0.05,0.05,0.001))
abline(v=0)
lines(PIX_DIFF$X50.,PIX_DIFF$ELEV,col="darkgreen",lwd=2)
polygon(x=c(PIX_DIFF$X25.,rev(PIX_DIFF$X75.)),y=c(PIX_DIFF$ELEV,rev(PIX_DIFF$ELEV)),
        col=adjustcolor("darkgreen",alpha.f = 0.5),border = "darkgreen")
