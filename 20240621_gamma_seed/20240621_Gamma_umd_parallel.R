
###############################################################################
#            Multivariate Density Estimation Using Bernstein Polynomials          #
#																	      	  #
#          Copyright(2021): Dan Han, University of Louisville, KY, USA
#																			  #	
#                         Version 2 - June 28, 2021	

#The following generates RMISE,ML1E and ML_infE for Gamma Distribution#	
###############################################################################




 setwd("C:/Users/Dan Han/OneDrive - University of Louisville/Research Shared with Collaborator/Multivariate Density Estimation Subject to Unimodality/20240621_gamma_seed")
# setwd("D:/OneDrive - University of Louisville/Research Shared with Collaborator/Multivariate Density Estimation Subject to Unimodality/20240621_gamma_seed")

library(mvtnorm)
library(MASS)
library(plyr)
library(lcmix)
library(sn)#generate skew-t distribution
library(ks)

library(kdevine)
source("Multivariate_denstiy_estimate.R")

########################################################
#The following generates RMISE for t Distribution#	
###########################################################

Sigma=matrix(c(4,6,6,16),ncol=2,nrow=2)
n=400
m=ceiling(n^0.4+1)
m=20
v0=100
alpha=c(2,2)



set.seed(300)


error_umd<- function (K){
  startTime<-Sys.time()
  RMISE_SUM=0
  ML1E_SUM=0
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


# Load necessary libraries
library(parallel)
library(doParallel)
library(foreach)

# Number of cores to use
numCores <- detectCores()-2


# Register the parallel backend
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Number of iterations
v0 <- 100 # Replace this with your actual v0 value

# Parallel loop
results <- foreach(v = 1:v0,.packages = c('quadprog', 'bivariate', 'magic', 'lcmix', 'MASS', 'plyr','sn', 'ks', 'kdevine')) %dopar% {
  cat("Starting iteration", v, "\n")
  umd_result <- error_umd(100)
  list(rmise = umd_result$RMISE_estimate,
       ml1e = umd_result$ML1E_estimate,
       ml_infe = umd_result$ML_infE_estimate)
}
abc=lapply(v=1:v0,FUN = error_umd(100))
# Ensure that results are in a list of lists
results <- lapply(results, function(res) {
  if (is.null(res)) {
    list(rmise = NA, ml1e = NA, ml_infe = NA)
  } else {
    res
  }
})

# Combine the results
rmise3_umd <- do.call(cbind, lapply(results, '[[', 'rmise'))
ML1E3_umd <- do.call(cbind, lapply(results, '[[', 'ml1e'))
ML_infE3_umd <- do.call(cbind, lapply(results, '[[', 'ml_infe'))

# Stop the cluster
stopCluster(cl)

# Unregister the parallel backend
registerDoSEQ()

# Print the combined results for verification
print(rmise3_umd)
print(ML1E3_umd)
print(ML_infE3_umd)




# rmise3_umd<-NULL
# ML1E3_umd<-NULL
# ML_infE3_umd<-NULL
# for(v in 1:v0){
#   umd_result=error_umd(100)
#   rmise3_umd=cbind(rmise3_umd,umd_result$RMISE_estimate)
#   ML1E3_umd=cbind(ML1E3_umd,umd_result$RMISE_estimate$ML1E_estimate)
#   ML_infE3_umd=cbind(ML_infE3_umd,umd_result$ML_infE_estimate)
# }

######################################################################

# load the data for the other methods, then draw the plot
df.rmise3_kde2d=data.frame(RMISE=t(rmise3_kde2d),group=rep("kde2d",v0))
df.rmise3_umd=data.frame(RMISE=t(rmise3_umd),group=rep("Berstein",v0))
df.rmise3_ks=data.frame(RMISE=t(rmise3_ks),group=rep("ks",v0))
df.rmise3_kdevine=data.frame(RMISE=t(rmise3_kdevine),group=rep("kdevine",v0))

rmise_gamma=rbind(df.rmise3_umd,df.rmise3_kde2d,df.rmise3_ks,df.rmise3_kdevine)

#draw boxplot for rmise t comparision
library(ggplot2)
boxplot_rmise_gamma=ggplot(rmise_gamma, aes(x=group,y=RMISE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(c(rmise3_umd,rmise3_kde2d,rmise3_ks,rmise3_kdevine))))

png(file = "boxplot_rmise_gamma.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_rmise_gamma
dev.off()

######################################################################

df.ML1E3_kde2d=data.frame(ML1E=t(ML1E3_kde2d),group=rep("kde2d",v0))
df.ML1E3_umd=data.frame(ML1E=t(ML1E3_umd),group=rep("Berstein",v0))
df.ML1E3_ks=data.frame(ML1E=t(ML1E3_ks),group=rep("ks",v0))
df.ML1E3_kdevine=data.frame(ML1E=t(ML1E3_kdevine),group=rep("kdevine",v0))

ML1E_gamma=rbind(df.ML1E3_umd,df.ML1E3_kde2d,df.ML1E3_ks,df.ML1E3_kdevine)

#draw boxplot for ML1E gamma comparision
library(ggplot2)
boxplot_ML1E_t=ggplot(ML1E_gamma, aes(x=group,y=ML1E)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML1E3_kde2d,ML1E3_umd,ML1E3_ks,ML1E3_kdevine)))

png(file = "boxplot_ML1E_gamma.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_ML1E_gamma
dev.off()

######################################################################


df.ML_infE3_kde2d=data.frame(ML_infE=t(ML_infE3_kde2d),group=rep("kde2d",v0))
df.ML_infE3_umd=data.frame(ML_infE=t(ML_infE3_umd),group=rep("Berstein",v0))
df.ML_infE3_ks=data.frame(ML_infE=t(ML_infE3_ks),group=rep("ks",v0))
df.ML_infE3_kdevine=data.frame(ML_infE=t(ML_infE3_kdevine),group=rep("kdevine",v0))

ML_infE_gamma=rbind(df.ML_infE3_umd,df.ML_infE3_kde2d,df.ML_infE3_ks,df.ML_infE3_kdevine)

#draw boxplot for ML_infE t comparision
library(ggplot2)
boxplot_ML_infE_gamma=ggplot(ML_infE_gamma, aes(x=group,y=ML_infE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML_infE3_kde2d, ML_infE3_umd,ML_infE3_ks,ML_infE3_kdevine)))

png(file = "boxplot_ML_inf_gamma.png", bg = "transparent",width=3000, height=2000, res=800, units="px")
boxplot_ML_infE_gamma
dev.off()

save.image(file="20240619_gamma_seed.RData")
