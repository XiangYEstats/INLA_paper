library(MASS)
library(sp)
library(gstat)
library(sf)
library(plyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(viridis)
library(INLA)
library(maptools)
library(automap)


# Data and source functions
load('aux/PollutantsMExicoCity2021.RData')
source('Functions.R')

# Mesh prep
CDMX_border<-st_union(CDMX_sf)
coords.sp.aux<-SpatialPoints(Stations[,5:6],proj4string = CRS("+proj=longlat +datum=WGS84"))
##transformed to KM
coords.sp<-spTransform(coords.sp.aux,CRS("+proj=utm +units=km +datum=WGS84 +zone=14"))
##Border
border.aux<- unionSpatialPolygons(CDMXpol,IDs = rep(1,16))
border<-spTransform(border.aux,CRS("+proj=utm +units=km +datum=WGS84 +zone=14"))

##Parameter space
# All parameters try all number in one decimal places.
# If possible, 2 decimal places would be better (but perhaps too intensive)
## max_edge_1: 5~20
## max_edge_2: 5~20
## min_angle_1: 10~30
## min_angle_2: 10~30
## cutoff: 1~10
## offset_1: 1~20
## offset_2: 1~20
## We start with three for now
me1<-seq(5,20,by=5)
me2<-seq(5,20,by=5)
co<-seq(1,10,by=2)

par.combined<-matrix(NA,nrow=length(me1)*length(me2)*length(co),ncol = 3)
aux.i<-1
for(i.me1 in 1:length(me1)){
  for(i.me2 in 1:length(me2)){
         for(i.co in 1:length(co)){
               par.combined[aux.i,]<-c(me1[i.me1],me2[i.me2],co[i.co])
              aux.i<-aux.i+1
            }}}
colnames(par.combined)<-c('edge1','edge2','cutoff')
n.pars<-nrow(par.combined)
##Testing data set
set.seed(89275)
M<-10
data.sim<-matrix(NA,nrow=M,ncol = length(coords.sp))
for(j in 1:M)data.sim[j,]<-GP.data(coords.sp,nu=1,kappa=4,sig=15,mu=45)

## now we run the Sim study
OUT.Sim<-vector('list',n.pars)
for( i.par in 1:n.pars){
  posteriors = as.data.frame(matrix(0,M,7))
  colnames(posteriors)=c("prior range", "prior sigma", "posterior range", "posterior sigma", "dic", "post_expect_se", "precision")
  for(j in 1:M){
    post = mesh.test2(data = data.sim[j,], coordinates=coords.sp,
                      boundary=border, max_edge_1= par.combined[i.par,1] , max_edge_2=par.combined[i.par,2],
                      min_angle_1=20, min_angle_2=20, cutoff=par.combined[i.par,3],
                      offset_1=5, offset_2=7)
    posteriors[j,1]=post$range0
    posteriors[j,2]=post$sigma0
    posteriors[j,3]=post$range_hat
    posteriors[j,4]=post$sigma_hat
    posteriors[j,5]=post$dic
    posteriors[j,6]=post$post_expect_se
    posteriors[j,7]=post$precision
  }
  OUT.Sim[[i.par]]<-posteriors
  print(i.par)
  }

save(file = 'aux/Output.RData',OUT.Sim,par.combined)
