
###############################################################################
#            Multivariate Density Estimation Using Bernstein Polynomials          #
#																	      	  #
#          Copyright(2021): Dan Han, University of Louisville, KY, USA
#																			  #	
#                         Version 2 - June 28, 2021	

#The following generates RMISE,ML1E and ML_infE for t Distribution#	
###############################################################################




setwd("C:/Users/danhan/OneDrive - University of Louisville/Research Shared with Collaborator/Multivariate Density Estimation Subject to Unimodality/20240619_t_seed")

library(mvtnorm)
library(MASS)
library(plyr)
library(lcmix)
library(sn)#generate skew-t distribution
# library(locpol)
library(ks)
# library(kdecopula)
# library(sparr)
library(kdevine)
source("Multivariate_denstiy_estimate_v2.R")
source("Multivariate_denstiy_estimate.R")

########################################################
#The following generates RMISE for t Distribution#	
###########################################################

Sigma=matrix(c(4,6,6,16),ncol=2,nrow=2)
n=400
m=ceiling(n^0.4+1)
m=30
v0=100
alpha=c(2,2)



set.seed(300)

error_umd<- function (K){
  startTime<-Sys.time()
  RMISE_SUM=0
  ML1E_SUM=0
  ML_infE_SUM=0
  for(k in 1:K){
    data=rmst(n=n, xi=rep(0,length(alpha)), Omega=Sigma, alpha=alpha,nu=6)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    result=umd_multi(data,crit="CN",m=m)
    density_estimate=result$dme_un(data)
    true_density=dmst(data,xi=rep(0,length(alpha)), Omega=Sigma, alpha=alpha,nu=6)
    diff_umd=density_estimate-true_density
    sum_e2=0
    sum_abe=0
    for(i in 1:n_dim){
      sum_e2=sum_e2+(abs(diff_umd[i]))^2
      sum_abe=sum_abe+abs(diff_umd[i])
      
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
    ML1E_SUM=ML1E_SUM+d*sum_abe
    
    sum_max_abs=max(abs(diff_umd))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
    
  }
  RMISE_estimate=RMISE_SUM/K
  ML1E_estimate=ML1E_SUM/K
  ML_infE_estimate=ML_infE_SUM/K
  endTime<-Sys.time()
  runningTime<-endTime-startTime
  return(list(RMISE_estimate=RMISE_estimate,ML1E_estimate=ML1E_estimate,ML_infE_estimate=ML_infE_estimate,runningTime=runningTime))}

rmise3_umd<-NULL
ML1E3_umd<-NULL
ML_infE3_umd<-NULL
for(v in 1:v0){
  umd_result=error_umd(100)
  rmise3_umd=cbind(rmise3_umd,umd_result$RMISE_estimate)
  ML1E3_umd=cbind(ML1E3_umd,umd_result$ML1E_estimate)
  ML_infE3_umd=cbind(ML_infE3_umd,umd_result$ML_infE_estimate)
}

######################################################################

# load the data for the other methods, then draw the plot
df.rmise3_kde2d=data.frame(RMISE=t(rmise3_kde2d),method=rep("kde2d",v0))
df.rmise3_umd=data.frame(RMISE=t(rmise3_umd),method=rep("Berstein",v0))
df.rmise3_ks=data.frame(RMISE=t(rmise3_ks),method=rep("ks",v0))
df.rmise3_kdevine=data.frame(RMISE=t(rmise3_kdevine),method=rep("kdevine",v0))

rmise_t=rbind(df.rmise3_umd,df.rmise3_kde2d,df.rmise3_ks,df.rmise3_kdevine)

#draw boxplot for rmise t comparision
library(ggplot2)
boxplot_rmise_t=ggplot(rmise_t, aes(x=method,y=RMISE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(c(rmise3_umd,rmise3_kde2d,rmise3_ks,rmise3_kdevine))))+
  theme(text = element_text(size = 15)) 

png(file = "boxplot_rmise_t_m30.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_rmise_t
dev.off()

######################################################################

df.ML1E3_kde2d=data.frame(ML1E=t(ML1E3_kde2d),method=rep("kde2d",v0))
df.ML1E3_umd=data.frame(ML1E=t(ML1E3_umd),method=rep("Berstein",v0))
df.ML1E3_ks=data.frame(ML1E=t(ML1E3_ks),method=rep("ks",v0))
df.ML1E3_kdevine=data.frame(ML1E=t(ML1E3_kdevine),method=rep("kdevine",v0))

ML1E_t=rbind(df.ML1E3_umd,df.ML1E3_kde2d,df.ML1E3_ks,df.ML1E3_kdevine)

#draw boxplot for ML1E t comparision
library(ggplot2)
boxplot_ML1E_t=ggplot(ML1E_t, aes(x=method,y=ML1E)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML1E3_kde2d,ML1E3_umd,ML1E3_ks,ML1E3_kdevine)))+
  theme(text = element_text(size = 15)) 

png(file = "boxplot_ML1E_t_m30.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_ML1E_t
dev.off()

######################################################################


df.ML_infE3_kde2d=data.frame(ML_infE=t(ML_infE3_kde2d),method=rep("kde2d",v0))
df.ML_infE3_umd=data.frame(ML_infE=t(ML_infE3_umd),method=rep("Berstein",v0))
df.ML_infE3_ks=data.frame(ML_infE=t(ML_infE3_ks),method=rep("ks",v0))
df.ML_infE3_kdevine=data.frame(ML_infE=t(ML_infE3_kdevine),method=rep("kdevine",v0))

ML_infE_t=rbind(df.ML_infE3_umd,df.ML_infE3_kde2d,df.ML_infE3_ks,df.ML_infE3_kdevine)

#draw boxplot for ML_infE t comparision
library(ggplot2)
boxplot_ML_infE_t=ggplot(ML_infE_t, aes(x=method,y=ML_infE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML_infE3_kde2d, ML_infE3_umd,ML_infE3_ks,ML_infE3_kdevine)))+
  theme(text = element_text(size = 15)) 

png(file = "boxplot_ML_inf_t_m30.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_ML_infE_t
dev.off()

save.image("20240619_t_seed_m30.RData")
