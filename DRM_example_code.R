library(BayesianTools)
library(mc2d)
library(truncnorm)
library(dplyr)
library(cubature)
library(raster)
library(rgdal)
library(sp) 
library(sf)

######################## LOAD IN ENVIRONMENTAL PREDICTORS #################################
Depth<-readRDS(file="Depth.rda")         #water depth             
SST-readRDS(file="SST.rda")              #sea surface temperature
maxHs<-readRDS(file="maxHs.rda")         # wave height

######################################## Dipseral Kernel function #######################################################
# Area to Area Dispersal function, Calculated based on Chipperfield et al. 2011 Methods in Ecology and Evolution 2011, 2, 668–681###
pDisp_LDD<-readRDS("pDisp_LDD.rda")

#########################  LOADING IN INITIAL SPECIES DENSITY PREDICTION ########################################
# Based on the results of correlatvie SDMs constructed from the initial year
N_matrix<-readRDS("N.rda")

#########################  LOADING IN SIZE CLASS OBSERVATIONS ########################################
y_SizeClass<-readRDS("y_SizeClass.rda")

#########################  LOADING IN PRESENCE-ABSENCE OBSERVATIONS ########################################
y_PresAbs<-readRDS("y_PresAbs.rda")


b<-0.001  # parameter for density dependence
Beta_p.1<-0
Beta_p.2<-28

#########################  DEFINE LIKELIHOOD FUNCTION ########################################


likelihood <- function (params) {
  Mort.1<-params[1]  
  Mort.2<-params[2]
  Mort.3<-params[3]
  Mort.4<-params[4]
  trans.G.1<-params[5]
  trans.G.2<-params[6]
  trans.S.1<-params[7]       
  trans.S.2<-params[8]  
  trans.F.3<-params[9]  
  trans.F.4<-params[10]  
  Beta_w.1<-params[11]
  Beta_w.2<-params[12]
  Beta_p.3<-params[13]
  Beta_w.3<-params[14]
  sigma.niche.1<-params[15]
  sigma.niche.2<-params[16]
  sigma.abundance<-params[17]
  
  
  for(k in 1:NumberOfYear)
  {
    
    ######################### Process 1: survival ########################################   
    u1[ ,1,k]<-rbinom(NoOfTotalGridCell,N_matrix[ ,1,k],(1-Mort.1))
    u1[ ,2,k]<-rbinom(NoOfTotalGridCell,N_matrix[ ,2,k],(1-Mort.2))
    u1[ ,3,k]<-rbinom(NoOfTotalGridCell,N_matrix[ ,3,k],(1-Mort.3))
    u1[ ,4,k]<-rbinom(NoOfTotalGridCell,N_matrix[ ,4,k],(1-Mort.4))
    
    
    ######################### Process 2: transition between size classes ##################
    ## Niche model - calculate transition rates 
    HS.mu[,k]<-exp(-((Depth-Beta_p.1)/Beta_w.1)^2-((SST[,k]-Beta_p.2)/Beta_w.2)^2)  
    HS[,k]<-rtruncnorm(NoOfTotalGridCell,a=0,b=1,mean=HS.mu[,k],sd=sigma.niche.1)
    f.g_[,k]<-HS[,k]*exp(-b*(N.total[,k]))  
    f.s[,k]<-1-HS[,k]

    
    ## transition between size classes in year k
    u2s1_a[ ,2,k]<-trans.G.1*f.g[ ,k]*u1[ ,1,k]        
    u2s1_a[ ,3,k]<-trans.G.2*f.g[ ,k]*u1[ ,1,k]       
    u2s1_a[ ,4,k]<-trans.G.2*f.g[ ,k]*u1[ ,1,k]
    u2s1.s<-u1[ ,1,k]-(u2s1_a[ ,2,k]+u2s1_a[ ,3,k]+u2s1_a[ ,4,k])
    u2s1_a[ ,1,k]<-ifelse(u2s1.s<=0,0.05,u2s1.s)   
    u2s1_p[ ,1:4,k]<-rdirichlet(NoOfTotalGridCell,u2s1_a[ ,1:4,k])   
    
    u2s2_a[ ,1,k]<-trans.S.1*f.s[ ,k]*u1[ ,2,k]        
    u2s2_a[ ,3,k]<-trans.G.1*f.g[ ,k]*u1[ ,2,k]       
    u2s2_a[ ,4,k]<-trans.G.2*f.g[ ,k]*u1[ ,2,k]
    u2s2.s<-u1_STC[ ,2,k]-(u2s2_a[ ,1,k]+u2s2_a[ ,3,k]+u2s2_a[ ,4,k])
    u2s2_a[ ,2,k]<-ifelse(u2s2.s<=0,0.05,u2s2.s)
    u2s2_p[ ,1:4,k]<-rdirichlet(NoOfTotalGridCell,u2s2_a[ ,1:4,k])  
    
    u2s3_a[ ,1,k]<-trans.S.2*f.s[ ,k]*u1[ ,3,k]      
    u2s3_a[ ,2,k]<-trans.S.1*f.s[ ,k]*u1[ ,3,k]     
    u2s3_a[ ,4,k]<-trans.G.1*f.g[ ,k]*u1[ ,3,k]
    u2s3.s<-u1[ ,3,k]-(u2s3_a[ ,1,k]+u2s3_a[ ,2,k]+u2s3_a[ ,4,k])
    u2s3_a[ ,3,k]<-ifelse(u2s3.s<=0,0.05,u2s3.s)
    u2s3_p[ ,1:4,k]<-rdirichlet(NoOfTotalGridCell,u2s3_a[ ,1:4,k])   
    
    u2s4_a[,1,k]<-trans.S.2*f.s[ ,k]*u1[ ,4,k]      
    u2s4_a[,2,k]<-trans.S.2*f.s[ ,k]*u1[ ,4,k]       
    u2s4_a[,3,k]<-trans.S.1*f.s[ ,k]*u1[ ,4,k] 
    u2s4.s<-u1_STC[ ,4,k]-(u2s4_a[ ,1,k]+u2s4_a[ ,2,k]+u2s4_a[ ,3,k])
    u2s4_a[,4,k]<-ifelse(u2s4.s<=0,0.05,u2s4.s)
    u2s4_p[,1:4,k]<-rdirichlet(NoOfTotalGridCell,u2s4_a[ ,1:4,k])   
    
    u2[,1,k]<-u1[,1,k]*u2s1_p[,1,k]+u1[,2,k]*u2s2_p[,1,k]+u1[,3,k]*u2s3_p[,1,k]+u1[,4,k]*u2s4_p[,1,k]
    u2[,2,k]<-u1[,1,k]*u2s1_p[,2,k]+u1[,2,k]*u2s2_p[,2,k]+u1[,3,k]*u2s3_p[,2,k]+u1[,4,k]*u2s4_p[,2,k]
    u2[,3,k]<-u1[,1,k]*u2s1_p[,3,k]+u1[,2,k]*u2s2_p[,3,k]+u1[,3,k]*u2s3_p[,3,k]+u1[,4,k]*u2s4_p[,3,k]
    u2[,4,k]<-u1[,1,k]*u2s1_p[,4,k]+u1[,2,k]*u2s2_p[,4,k]+u1[,3,k]*u2s3_p[,4,k]+u1[,4,k]*u2s4_p[,4,k]
    
    
    ######################### Process 3: fragmentation ###########################
    #calculate how many fragments generated in each cell 
    
    f.f.mu[,k]<-exp(-((maxHs[,k]-Beta_p.3)/Beta_w.3)^2)
    f.f[,k]<-rtruncnorm(NoOfTotalGridCell,a=0,b=1,mean=f.f.mu[,k],sd=sigma.niche.2)
    u3_frag.vector[,k]<-f.f[,k]*trans.F.3*u2[,3,k]+f.f[,k]*trans.F.4*u2[,4,k] 
  
    ######################### Process 4: fragments dispersal ###########################    
    ## Convert vector to raster for dispersal moving winodow calculation
    u3_frag.raster<-raster(nrow=NoOfRows, ncol=NoOfColumna,ext=extent(maxHs),crs=crs(maxHs),resolution=res(maxHs))
    u3_frag.raster[]<-u3_frag.vector[,k]
    u3_frag.raster<-raster::mask(u3_frag.raster,maxHsk)
    
    ## Moving windows
    u3_frag.raster.mv<-focal(u3_frag.raster,w=pDisp_LDD,fun=sum,na.rm=T)
    
    ## Convert back to vectors
    # dispersal probability is the same for size class 1 and 2 fragments
    u4_frag_post[,1,k]<-as.vector(u3_frag.raster.mv) 
    u4_frag_post[,2,k]<-as.vector(u3_frag.raster.mv) 
    
    #Convert all the NAs (masked out grids) to 0 to avoid later summation problem
    u4_frag_post[,1,k]<-ifelse(is.na(u4_frag_post[,1,k])==F,u4_frag_post[,1,k],0) 
    u4_frag_post[,2,k]<-ifelse(is.na(u4_frag_post[,2,k])==F,u4_frag_post[,2,k],0) 
    
    
    ######################### Process 5: final summation ###########################################
    N[ ,1,k+1]<-round(u2[,1,k]+u4_frag_post[,1,k])
    N[ ,2,k+1]<-round(u2[,2,k]+u4_frag_post[,2,k])
    u4_frag_post_s3<-ifelse(is.na(f.f[,k]*trans.F.4*u2[,4,k])==F,(f.f[,k]*trans.F.4*u2[,4,k]),0)
    N[ ,3,k+1]<-round(u2[,3,k]+u4_frag_post_s3)
    N[ ,4,k+1]<-round(u2[,4,k])
    
    N.total[,k+1]<-N[ ,1,k+1]+N[ ,2,k+1]+N[ ,3,k+1]+ N[ ,4,k+1]
  }
  
  
  LogN<-log(N+0.00001)
  LogN.total<-log(N.total+0.00001)
  
  ######################### Likelihood/Data model ###########################################
  
  # Data source 1: composition data
  
  for (i in 1:length(SizeClassDataGridNumber)) {
    
    lambda.1<-rlnorm(1,meanlog=LogN[DataGridNumber[i],1],sdlog=sigma.abundance)
    lambda.2<-rlnorm(1,meanlog=LogN[DataGridNumber[i],2],sdlog=sigma.abundance)
    lambda.3<-rlnorm(1,meanlog=LogN[DataGridNumber[i],3],sdlog=sigma.abundance)
    lambda.4<-rlnorm(1,meanlog=LogN[DataGridNumber[i],4],sdlog=sigma.abundance)
    lambda.sum<-lambda.1+lambda.2+lambda.3+lambda.4
    
    if (lambda.sum==0&&y_SizeClass[i,6]>0) {
      LikSum[i]<-dpois(y_SizeClass[i,6],lambda.sum ,log=T)
    }
    
    if (lambda.sum>0&&y_SizeClass[i,6]>0) {
      comp.p[i,1]<-lambda.1/lambda.sum 
      comp.p[i,2]<-lambda.2/lambda.sum
      comp.p[i,3]<-lambda.3/lambda.sum
      comp.p[i,4]<-lambda.4/lambda.sum
      
      LikSum[i]<-dmultinom(y_SizeClass[i,2:5],size=y_SizeClass[i,6],prob=comp.p[i,1:4],log=TRUE)     
    }
   }
  
  # Data source 2: presence/absence data
  
  for(i in 1:length(PresAbsDataGridNumber)){
    pi.abs<-rbeta(1,90,10)  #Obtained from Pagel et al 2012 
    phi<-1-pi.abs^N.total[PresAbsDataGridNumber[i],2]+10^-5  
    LikPres[i]<-dbinom(y_PresAbs[i,2], size=1,prob=phi, log = T) 
  }
  
  LL<-sum(LikSum,LikPres)
  return(LL)
}
##############################################################


########## Define prior function ##############################
PriorDensity<-function(params){
  Mort.1<-params[1]  
  Mort.2<-params[2]
  Mort.3<-params[3]
  Mort.4<-params[4]
  trans.G.1<-params[5]
  trans.G.2<-params[6]
  trans.S.1<-params[7]       
  trans.S.2<-params[8]  
  trans.F.3<-params[9]  
  trans.F.4<-params[10]  
  Beta_w.1<-params[11]
  Beta_w.2<-params[12]
  Beta_p.3<-params[13]
  Beta_w.3<-params[14]
  sigma.niche.1<-params[15]
  sigma.niche.2<-params[16]
  sigma.abundance<-params[17]
  
  Mort.1.prior<-dbeta(Mort.1,0.0129338, 1.162866,log=T)
  Mort.2.prior<-dbeta(Mort.2,2.2287, 10.8813,log=T)
  Mort.3.prior<-dbeta(Mort.3,0.2400384,4.046362,log=T)
  Mort.4.prior<-dbeta(Mort.4,0.0129338, 1.162866,log=T)
  trans.G.1.prior<-dbeta(trans.G.1,1.688948,9.826607,log=T)      
  trans.G.2.prior<-dbeta(trans.G.2,0.0098,0.9702,log=T) 
  trans.S.1.prior<-dbeta(trans.S.1,0.6471,6.5429,log=T)
  trans.S.2.prior<-dbeta(trans.S.2,0.0098,0.9702,log=T)
  trans.F.3.prior<-dbeta(trans.F.3,0.1875,3.5625,log=T)
  trans.F.4.prior<-dbeta(trans.F.4,3.320286,12.3661,log=T)
  Beta_w.1.prior<-dgamma(Beta_w.1,900, 60,log=T)    
  Beta_w.2.prior<-log(dtruncnorm(Beta_w.2,a=10^-5,b=100, mean=1,sd=1))           
  Beta_p.3.prior<-log(dtruncnorm(Beta_p.3,a=10^-5,b=100, mean=1.206228,sd=0.5))   
  Beta_w.3.prior<-log(dtruncnorm(Beta_w.3,a=10^-5,b=100, mean=1,sd=1))    
  sigma.niche.1.prior<-dgamma(sigma.niche.1,0.01,0.01,log=T)
  sigma.niche.2.prior<-dgamma(sigma.niche.2,0.01,0.01,log=T)
  sigma.abundance.prior<-dgamma(sigma.abundance,0.01,0.01,log=T)
  
  prior.sum<- sum(Mort.1.prior,Mort.2.prior,Mort.3.prior,Mort.4.prior,trans.G.1.prior,trans.G.2.prior,trans.S.1.prior,
                  trans.S.2.prior,trans.F.3.prior,trans.F.4.prior,Beta_w.1.prior,
                  Beta_w.2.prior,Beta_p.3.prior,Beta_w.3.prior,sigma.niche.1.prior,sigma.niche.2.prior,sigma.abundance.prior)  
  return(prior.sum) 
}

PriorSampler<-function(n=1){
  Mort.1<-rbeta(n,0.0129338, 1.162866)
  Mort.2<-rbeta(n,2.2287, 10.8813)
  Mort.3<-rbeta(n,0.2400384,4.046362)
  Mort.4<-rbeta(n,0.0129338, 1.162866)
  trans.G.1<-rbeta(n,1.688948,9.826607)    
  trans.G.2<-rbeta(n,0.0098,0.9702)  
  trans.S.1<-rbeta(n,0.6471,6.5429)
  trans.S.2<-rbeta(n,0.0098,0.9702)
  trans.F.3<-rbeta(n,0.1875,3.5625)
  trans.F.4<-rbeta(n,3.320286,12.3661)
  Beta_w.1<-rgamma(n,900, 60)   
  Beta_w.2<-rtruncnorm(n,a=10^-5,b=100, mean=1,sd=1)         
  Beta_p.3<-rtruncnorm(n,a=10^-5,b=100, mean=1.206228,sd=0.5)
  Beta_w.3<-rtruncnorm(n,a=10^-5,b=100, mean=1,sd=1)     
  sigma.niche.1<-rgamma(n,0.01,0.01)
  sigma.niche.2<-rgamma(n,0.01,0.01)
  sigma.abundance<-rgamma(n,0.01,0.01)
  return(cbind(Mort.1,Mort.2,Mort.3,Mort.4,trans.G.1,trans.G.2,trans.S.1,trans.S.2,trans.F.3,trans.F.4,Beta_w.1,Beta_w.2,Beta_p.3,Beta_w.3,sigma.niche.1,sigma.niche.2,sigma.abundance)) 
}


Prior<- createPrior(density = PriorDensity, sampler = PriorSampler,lower=c(rep(0.00001,17)),upper=c(rep(35,17)))
bayesianSetup <- createBayesianSetup(likelihood=likelihood,prior=Prior,parallel=10)
settings <- list(iterations =NumberOfIterations,consoleUpdates=1,nrChains = 2,startValue =initialValues,burnin=NumberOfBurnIn)
Results.Coral.DRM<- runMCMC(bayesianSetup = bayesianSetup, sampler = "DEzs", settings = settings)
