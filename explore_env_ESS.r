

# install.packages("SpatialPack",dep=T)

library(SpatialPack)
library(dplyr)
library(readr)
library(ape)
library(spdep)
library(leaflet)

setwd("~/STUDENTS and COLLABORATIONS/Eucalyptus GEA workshop/euclandscape-analysisgroup")

datadir="rawdata"
envfilename="JordanEmicroMetadataBySample.csv"
envfile=file.path(datadir,envfilename)

## read data
envsamp=read_csv(envfile)
envnames_all=names(envsamp[-c(1:2)])
env=envsamp %>% 
  group_by(Population) %>% 
  summarise_at(.funs=list(mean),.vars=vars(envnames_all))
env[,envnames_all]=scale(env[envnames_all])

#########################################################################
## Pairwise analyses
#########################################################################
## Define response and predictor
envnames=envnames_all[which(!envnames_all %in% c("Lat","Lon"))]

# extracting the coordinates from the radiata data set
coords <- as.data.frame(env[c("Lon","Lat")])

## set up output data frame
mymatrix=matrix(nrow=length(envnames),ncol=length(envnames)) 
rownames(mymatrix)=envnames
colnames(mymatrix)=envnames
mymatrix2=mymatrix

for(i in 1:length(envnames)){
  # i=2
  xname=envnames[i]
  x=t(t(env[,xname]))
  for(j in c(1:length(envnames))[which(!envnames==xname)]){
    yname=envnames[j]
    y=unlist(env[yname])
    
    # computing the modified F-test of spatial association
    z <- modified.Ftest(x, y, coords)
    mymatrix[i,j]=z$ESS
    
    # computing the modified t-test of spatial association
    z <- modified.ttest(x, y, coords)
    mymatrix2[i,j]=z$ESS
    
  }
}
mymatrix
mymatrix2

## ESS is the same for the F-test with two variables as for the t-test of a correlation

#########################################################################
## Autocorrelation - Moran's I in ape and APLE in spdep
#########################################################################

## Set up frame to hold output
envI=data.frame(env=envnames,MoransI=NA,pval=NA,APLE=NA)

## weights for ape
w=as.matrix(1/dist(coords,upper = T))
## weights for spdep
listw0=mat2listw(as.matrix(1/dist(coords)))
listw1=nb2listw(listw0$neighbours,listw0$weights,"W")

for(i in 1:length(envnames)){
  varname=envnames[i]
  x=unlist(env[,varname])
  # x=c(scale(unlist(env[,varname])))
  
  y=Moran.I(x,w)
  envI$MoransI[i]=y$observed
  envI$pval[i]=y$p.value
  envI$APLE[i]=aple(x,listw1)
  
}

envI


#########################################################################
## Effective sample size for the mean (Griffiths 2005)
#########################################################################

## n is the sample size
## rho is the spatial autocorrelation coefficient for a SAR model specification
##    - using Moran's I, as APLE gives negative values
ess=function(n,rho){
  n_ess=n*(1-(1/(1-exp(-1.92369))*((n-1)/n)*(1-exp(-2.12373*rho+0.20024*sqrt(rho)))))
  n_ess
}

ess(26,envI$MoransI)

