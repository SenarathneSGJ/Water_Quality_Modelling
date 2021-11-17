library(riverdist)
library(rgdal)
library(raster)
library(mvtnorm)
library(fields)
library(corpcor)
library(Matrix)

source("Univar_loglike.R")
source("Univar_pred.R")
source("Tivar_loglike.R")
source("Tivar_pred.R")

### load the dataset which contain monthly means
Data.ann = read.csv("RiverData_Monthly.csv",header = TRUE)
Data.ann = Data.ann[order(Data.ann$Year,Data.ann$Month), ]
#View(Data.ann)

Time= (Data.ann$Year-2010)*12+Data.ann$Month
Data.ann=cbind(Data.ann,Time)
Data.ann=Data.ann[Data.ann$Year!=2010,]

Data1= Data.ann[,c(2,19:20)]
Data2=unique(Data1)

#### Creare the distance matrix
source("Custom_riverdist_func.R")
MySHP = "MySHP/HydroRIVERS_v10_eu.shp"

AOI = readOGR(MySHP)
crs(AOI)

Proj.shape = "+proj=utm +zone=30 +datum=WGS84 +no_defs +ellps=WGS84"
MyRivernetwork <- my.custom.line2network(path=MySHP,layer="HydroRIVERS_v10_eu",tolerance = 100,reproject=Proj.shape)

cord.dec = SpatialPoints(cbind(Data2$Longitude, Data2$Latitude), proj4string=CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS(Proj.shape))
locs=as.data.frame(cord.UTM)

#### convert XY (stations) to river location
River.XY <- xy2segvert(x=locs$coords.x1, y=locs$coords.x2, rivers=MyRivernetwork)

#Computing distance matrix between all observations
dmat <- riverdistancemat(River.XY$seg,River.XY$vert,MyRivernetwork)
dist.Y <- dmat/max(dmat)

site_name=unique(Data.ann$Site)
site_no=1:length(site_name)
d.site=data.frame(site_name,site_no)

### Start

for(j in 1:500){

set.seed(1231*j)

pred.site=d.site[sort(sample(nrow(d.site),5,replace=F)),] # For n=5 prediction locations

samp_last=Data.ann[441:462,]
samp_last2=samp_last[!samp_last$Site%in%pred.site$site_name,]
samp_no=which(!samp_last$Site%in%pred.site$site_name)
Pred_loc_data=samp_last[samp_last$Site%in%pred.site$site_name,]

Data.ann2= rbind(Data.ann[1:440,],samp_last2)

r01 <- 0.001
r02 <- 0.001
r03 <- 0.001

dataX=(Data.ann2[,c(7,8,21)])
dataY=Data.ann2[,c(6,17,14)]

num.knots=3
TimeD=c(scale(dataX[,3]))

knots <- quantile(unique(TimeD),
                  seq(0, 1, length = (num.knots + 2))[-c(1, (num.knots + 2))])
Z_K <- (abs(outer(TimeD, knots, "-"))) ^ 3
OMEGA_all <- (abs(outer(knots, knots, "-"))) ^ 3
svd.OMEGA_all <- svd(OMEGA_all)
sqrt.OMEGA_all <- t(svd.OMEGA_all$v %*%
                      (t(svd.OMEGA_all$u) * sqrt(svd.OMEGA_all$d)))
Z <- t(solve(sqrt.OMEGA_all, t(Z_K)))

#Estimate the parameters of the model using MLE
Theta.init <- c(40,-13,0.1,-3,10,-5,-0.1,-0.1,40,-10,0.1,0.1,1,1,1,1,1,1,5,5,5,0,0,0,.1,.1,.1,.1,.1,.1,.1,.1,.1) #Initial values for nlm function (trivariate model)

lp.approx <- nlm(p=Theta.init,samp_no=samp_no,dataX=log(Data.ann2[,c(7,8,21)]),dist.Y=dist.Y,dataY=Data.ann2[,c(6,17,14)],z=Z, f=NLog_Like, hessian=TRUE,iterlim=400) 
Post_mean=lp.approx$estimate

Theta.init1 <- c(40,-13,0.01,-3,1,1,3,-.1,-.1,-.1) #initial values for nlm function (univariate model with Y1)
Theta.init2 <- c(10,-10,0.1,0.1,1,1,3,.1,.1,.1) #initial values for nlm function (univariate model with Y2)
Theta.init3 <- c(10,-10,0.1,0.1,1,1,3,.1,.1,.1) #initial values for nlm function (univariate model with Y3)

lp.approx2 <- nlm(p=Theta.init1,samp_no=samp_no,dataX=log(Data.ann2[,c(7,8,21)]),dist.Y=dist.Y,dataY=data.frame(Data.ann2[,c(6)]),z=Z, f=NLog_Like_Res1, hessian=TRUE,iterlim=400)
Post_mean1=lp.approx2$estimate

lp.approx3 <- nlm(p=Theta.init2,samp_no=samp_no,dataX=log(Data.ann2[,c(7,8,21)]),dist.Y=dist.Y,dataY=data.frame(Data.ann2[,c(17)]),z=Z, f=NLog_Like_Res1, hessian=TRUE,iterlim=400)
Post_mean2=lp.approx3$estimate

lp.approx4 <- nlm(p=Theta.init3,samp_no=samp_no,dataX=log(Data.ann2[,c(7,8,21)]),dist.Y=dist.Y,dataY=data.frame(Data.ann2[,c(14)]),z=Z, f=NLog_Like_Res1, hessian=TRUE,iterlim=400)
Post_mean3=lp.approx4$estimate

#### predictions 
TR.para = Post_mean
TR.para1 = Post_mean1
TR.para2 = Post_mean2
TR.para3 = Post_mean3
Fac.loc= Pred_loc_data[,c(7,8,21)]
loc.No=which(samp_last$Site%in%Pred_loc_data$Site)
dist.mat= dmat[loc.No,loc.No]/max(dmat)

z_new=tail(Z,5) # For n=5 prediction locations

pred.joint=matrix(ncol=3,nrow=100)
predY1=matrix(ncol=nrow(Fac.loc),nrow=100)
predY2=matrix(ncol=nrow(Fac.loc),nrow=100)
predY3=matrix(ncol=nrow(Fac.loc),nrow=100)
predY1U=matrix(ncol=nrow(Fac.loc),nrow=100)
predY2U=matrix(ncol=nrow(Fac.loc),nrow=100)
predY3U=matrix(ncol=nrow(Fac.loc),nrow=100)

for(i in 1:100){

Pred.Y=pred(TR.para,Fac.loc=log(Fac.loc),dist.mat,z_new=z_new)
Pred.Y1=pred2(TR.para1,Fac.loc=log(Fac.loc),dist.mat,z_new=z_new)
Pred.Y2=pred2(TR.para2,Fac.loc=log(Fac.loc),dist.mat,z_new=z_new)
Pred.Y3=pred2(TR.para3,Fac.loc=log(Fac.loc),dist.mat,z_new=z_new)

predY1[i,]= Pred.Y[,1]
predY2[i,]= Pred.Y[,2]
predY3[i,]= Pred.Y[,3]

predY1U[i,]= Pred.Y1[,1]
predY2U[i,]= Pred.Y2[,1]
predY3U[i,]= Pred.Y3[,1]

Act.Y=Pred_loc_data[,c(6,17,14)]
#plot(Pred.Y[,1],Act.Y[,1])
cor(Pred.Y[,1],Act.Y[,1])

#plot(Pred.Y[,2],Act.Y[,2])
cor(Pred.Y[,2],Act.Y[,2])

pred.error=colMeans((Act.Y-Pred.Y)^2)
pred.joint[i,]=pred.error
}

MY1=apply(predY1,2,mean)
MY2=apply(predY2,2,mean)
MY3=apply(predY3,2,mean)

MY1U=apply(predY1U,2,mean)
MY2U=apply(predY2U,2,mean)
MY3U=apply(predY3U,2,mean)

SY1=apply(predY1,2,sd)
SY2=apply(predY2,2,sd)
SY3=apply(predY3,2,sd)

SY1U=apply(predY1U,2,sd)
SY2U=apply(predY2U,2,sd)
SY3U=apply(predY3U,2,sd)

data1=data.frame(Station=as.factor(pred.site$site_no),X=pred.site$site_name,MY1,SY1,MY1U,SY1U,MY2U,SY2U,MY2,SY2,MY3,SY3,MY3U,SY3U,Y1=Act.Y[,1],Y2=Act.Y[,2],Y3=Act.Y[,3])

vec.name=paste("Results2/out",j,".RData",sep="")
save(data1,file=vec.name)
}
