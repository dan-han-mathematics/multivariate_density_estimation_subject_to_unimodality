
###############################################################################
#            Multivariate Density Estimation Using Bernstein Polynomials          #
#																	      	  #
#          Copyright(2021): Dan Han, University of Louisville, KY, USA
#																			  #	
#                         Version 2 - June 28, 2021	

#The following generates RMISE,ML1E and ML_infE for Gamma Distribution#	
###############################################################################




setwd("C:/Users/Dan Han/OneDrive - University of Louisville/Research Shared with Collaborator/Multivariate Density Estimation Subject to Unimodality/20240528_Gamma_seed")

library(mvtnorm)
library(MASS)
library(plyr)
library(lcmix)
library(ggplot2)
library(dplyr)
library(ks)
library(kdevine)
source("Multivariate_denstiy_estimate.R")

########################################################
#The following generates RMISE for Gamma Distribution#	
###########################################################

Sigma=matrix(c(4,6,6,16),ncol=2,nrow=2)
n=400
m=ceiling(n^0.4+1)
m=10
v0=100

set.seed(200)

RMISE_umd<- function (K){
  RMISE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    result=umd_multi(data,crit="CN",m=m)
    density_estimate=result$dme_un(data)
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_umd=density_estimate-true_density
    sum_e2=0
    for(i in 1:n_dim){
      sum_e2=sum_e2+(abs(diff_umd[i]))^2
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
  }
  RMISE_estimate=RMISE_SUM/K
  return(RMISE_estimate)}


RMISE_kde2d<- function (K){
  RMISE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    mass_result=kde2d(x=data[,1],y=data[,2],n = 20, lims = c(range(data[,1]), range(data[,2])))
    
    # create a new data frame of that 2d density grid
    gr <- data.frame(with(mass_result, expand.grid(x,y)), as.vector(mass_result$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    
    # Apply the model to the original data to estimate density at that point
    smoothed_density<- predict(mod, newdata=data.frame(xgr=data[,1], ygr=data[,2]))
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_kde2d=as.vector(smoothed_density-true_density)
    sum_e2=0
    for(i in 1:n_dim){
      sum_e2=sum_e2+(abs(diff_kde2d[i]))^2
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
  }
  RMISE_estimate=RMISE_SUM/K
  return(RMISE_estimate)}


RMISE_ks<- function (K){
  RMISE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    Hpi1 <- Hpi(x=data)
    ks_result= ks::kde(x=data, H=Hpi1,eval.points = data)
    density_estimate=ks_result$estimate
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_ks=density_estimate-true_density
    sum_e2=0
    for(i in 1:n_dim){
      sum_e2=sum_e2+(abs(diff_ks[i]))^2
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
  }
  RMISE_estimate=RMISE_SUM/K
  return(RMISE_estimate)}

RMISE_kdevine<- function (K){
  RMISE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    fit <- kdevine(data)
    density_estimate=dkdevine(data, fit)
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_kdevine=density_estimate-true_density
    sum_e2=0
    for(i in 1:n_dim){
      sum_e2=sum_e2+(abs(diff_kdevine[i]))^2
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
  }
  RMISE_estimate=RMISE_SUM/K
  return(RMISE_estimate)}


rmise3_umd<-NULL
rmise3_kde2d<-NULL
rmise3_ks<-NULL
rmise3_kdevine<-NULL
for(v in 1:v0){
  rmise3_umd=cbind(rmise3_umd,RMISE_umd(100))
  rmise3_kde2d=cbind(rmise3_kde2d,RMISE_kde2d(100))
  rmise3_ks=cbind(rmise3_ks,RMISE_ks(100))
  rmise3_kdevine=cbind(rmise3_kdevine,RMISE_kdevine(100))
  
}

df.rmise3_kde2d=data.frame(RMISE=t(rmise3_kde2d),group=rep("kde2d",v0))
df.rmise3_umd=data.frame(RMISE=t(rmise3_umd),group=rep("Berstein",v0))
df.rmise3_ks=data.frame(RMISE=t(rmise3_ks),group=rep("ks",v0))
df.rmise3_kdevine=data.frame(RMISE=t(rmise3_kdevine),group=rep("kdevine",v0))

rmise_gamma=rbind(df.rmise3_umd,df.rmise3_kde2d,df.rmise3_ks,df.rmise3_kdevine)



#draw boxplot for rmise gamma comparision
library(ggplot2)
boxplot_rmise_gamma=ggplot(rmise_gamma, aes(x=group,y=RMISE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(c(rmise3_umd,rmise3_kde2d,rmise3_ks,rmise3_kdevine))))


png(file = "boxplot_rmise_gamma.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_rmise_gamma
dev.off()

boxplot(t(rmise3_umd))
boxplot(t(rmise3_kde2d))


########################################################
#The following generates ML1E for Gamma Distribution#	
###########################################################
ML1E_umd<- function (K){
  ML1E_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    result=umd_multi(data,crit="CN",m=m)
    density_estimate=result$dme_un(data)
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_umd=density_estimate-true_density
    sum_abe=0
    for(i in 1:n_dim){
      sum_abe=sum_abe+abs(diff_umd[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}


ML1E_kde2d<- function (K){
  ML1E_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    mass_result=kde2d(x=data[,1],y=data[,2],n = 20, lims = c(range(data[,1]), range(data[,2])))
    
    # create a new data frame of that 2d density grid
    gr <- data.frame(with(mass_result, expand.grid(x,y)), as.vector(mass_result$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    
    # Apply the model to the original data to estimate density at that point
    smoothed_density<- predict(mod, newdata=data.frame(xgr=data[,1], ygr=data[,2]))
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_kde2d=as.vector(smoothed_density-true_density)
    sum_abe=0
    for(i in 1:n_dim){
      sum_abe=sum_abe+abs(diff_kde2d[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}

ML1E_ks<- function (K){
  ML1E_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    Hpi1 <- Hpi(x=data)
    ks_result= ks::kde(x=data, H=Hpi1,eval.points = data)
    density_estimate=ks_result$estimate
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_ks=density_estimate-true_density
    sum_abe=0
    for(i in 1:n_dim){
      sum_abe=sum_abe+abs(diff_ks[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}

ML1E_kdevine<- function (K){
  ML1E_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    fit <- kdevine(data)
    density_estimate=dkdevine(data, fit)
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_kdevine=density_estimate-true_density
    sum_abe=0
    for(i in 1:n_dim){
      sum_abe=sum_abe+abs(diff_kdevine[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}

ML1E3_umd<-NULL
ML1E3_kde2d<-NULL
ML1E3_ks<-NULL
ML1E3_kdevine<-NULL
for(v in 1:v0){
  ML1E3_umd=cbind(ML1E3_umd,ML1E_umd(100))
  ML1E3_kde2d=cbind(ML1E3_kde2d,ML1E_kde2d(100))
  ML1E3_ks=cbind(ML1E3_ks,ML1E_ks(100))
  ML1E3_kdevine=cbind(ML1E3_kdevine,ML1E_kdevine(100))
  
  
}



df.ML1E3_kde2d=data.frame(ML1E=t(ML1E3_kde2d),group=rep("kde2d",v0))
df.ML1E3_umd=data.frame(ML1E=t(ML1E3_umd),group=rep("Berstein",v0))
df.ML1E3_ks=data.frame(ML1E=t(ML1E3_ks),group=rep("ks",v0))
df.ML1E3_kdevine=data.frame(ML1E=t(ML1E3_kdevine),group=rep("kdevine",v0))

ML1E_gamma=rbind(df.ML1E3_umd,df.ML1E3_kde2d,df.ML1E3_ks,df.ML1E3_kdevine)

#draw boxplot for ML1E Gamma comparision
library(ggplot2)
boxplot_ML1E_gamma=ggplot(ML1E_gamma, aes(x=group,y=ML1E)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML1E3_kde2d,ML1E3_ks,ML1E3_kdevine)))


png(file = "boxplot_ML1E_gamma.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_ML1E_gamma
dev.off()

boxplot(t(ML1E3_umd))
boxplot(t(ML1E3_kde2d))

########################################################
#The following generates ML_infE for Gamma Distribution#	
###########################################################
ML_infE_umd<- function (K){
  ML_infE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    result=umd_multi(data,crit="CN",m=m)
    density_estimate=result$dme_un(data)
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_umd=density_estimate-true_density
    sum_max_abs=max(abs(diff_umd))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}


ML_infE_kde2d<- function (K){
  ML_infE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    
    mass_result=kde2d(x=data[,1],y=data[,2],n = 20, lims = c(range(data[,1]), range(data[,2])))

    # create a new data frame of that 2d density grid
    gr <- data.frame(with(mass_result, expand.grid(x,y)), as.vector(mass_result$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    
    # Apply the model to the original data to estimate density at that point
    smoothed_density<- predict(mod, newdata=data.frame(xgr=data[,1], ygr=data[,2]))
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_kde2d=as.vector(smoothed_density-true_density)
    sum_max_abs=max(abs(diff_kde2d))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}



ML_infE_ks<- function (K){
  ML_infE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    Hpi1 <- Hpi(x=data)
    ks_result=ks::kde(x=data, H=Hpi1,eval.points = data)
    density_estimate=ks_result$estimate
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_ks=density_estimate-true_density
    sum_max_abs=max(abs(diff_ks))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}

ML_infE_kdevine<- function (K){
  ML_infE_SUM=0
  for(k in 1:K){
    data=rmvgamma(n=n,shape=c(4,6),rate=c(1,1),corr=Sigma)
    data=data[is.finite(data[,2]) == 1 & is.finite(data[,1])==1, ]
    
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    fit <- kdevine(data)
    density_estimate=dkdevine(data, fit)
    true_density=dmvgamma(data,shape=c(4,6),rate=c(1,1),corr=Sigma,log=FALSE)
    diff_kdevine=density_estimate-true_density
    sum_max_abs=max(abs(diff_kdevine))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}

ML_infE3_umd<-NULL
ML_infE3_kde2d<-NULL
ML_infE3_ks<-NULL
ML_infE3_kdevine<-NULL
for(v in 1:v0){
  ML_infE3_umd=cbind(ML_infE3_umd,ML_infE_umd(100))
  ML_infE3_kde2d=cbind(ML_infE3_kde2d,ML_infE_kde2d(100))
  ML_infE3_ks=cbind(ML_infE3_ks,ML_infE_ks(100))
  ML_infE3_kdevine=cbind(ML_infE3_kdevine,ML_infE_kdevine(100))
  
}

df.ML_infE3_kde2d=data.frame(ML_infE=t(ML_infE3_kde2d),group=rep("kde2d",v0))
df.ML_infE3_umd=data.frame(ML_infE=t(ML_infE3_umd),group=rep("Berstein",v0))
df.ML_infE3_ks=data.frame(ML_infE=t(ML_infE3_ks),group=rep("ks",v0))
df.ML_infE3_kdevine=data.frame(ML_infE=t(ML_infE3_kdevine),group=rep("kdevine",v0))

ML_infE_gamma=rbind(df.ML_infE3_umd,df.ML_infE3_kde2d,df.ML_infE3_ks,df.ML_infE3_kdevine)


#draw boxplot for ML_infE Gamma comparision
library(ggplot2)
boxplot_ML_infE_gamma=ggplot(ML_infE_gamma, aes(x=group,y=ML_infE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML_infE3_kde2d,ML_infE3_ks,ML_infE3_kdevine)))


png(file = "boxplot_ML_inf_gamma.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_ML_infE_gamma
dev.off()

boxplot(t(ML_infE3_umd))
boxplot(t(ML_infE3_kde2d))

save.image("20240528_Gamma_seed_result.RData")
