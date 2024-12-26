
###############################################################################
#            Multivariate Density Estimation Using Bernstein Polynomials          #
#																	      	  #
#          Copyright(2021): Dan Han, University of Louisville, KY, USA
#																			  #	
#                         Version 2 - June 28, 2021	

#The following generates RMISE,ML1E and ML_infE for Gaussian Distribution#	
###############################################################################



#dell laptop   
# setwd("C:/Users/DanHan/OneDrive - University of Louisville/Research Shared with Collaborator/Multivariate Density Estimation Subject to Unimodality/20240610_normal_outlier")
#surface laptop
setwd("C:/Users/Dan Han/OneDrive - University of Louisville/Research Shared with Collaborator/Multivariate Density Estimation Subject to Unimodality/20240610_normal_outlier")

library(mvtnorm)
library(MASS)
# library(raster)
library(plyr)
library(ks)
library("Rlab")    
library(kdevine)
library(ggplot2)
library(dplyr)
library(reshape2)

source("Multivariate_denstiy_estimate.R")
########################################################
#The following generates RMISE for Gaussian Distribution#	
###########################################################

Sigma=matrix(c(4,6,6,16),ncol=2,nrow=2)
n=400
m=ceiling(n^0.4+1)
m=10
v0=100

set.seed(100)


RMISE_umd<- function (K,data){
        RMISE_SUM=0
        for(k in 1:K){
                d=((max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2])))/(n^2)
                result=umd_multi(data,crit="CN",m=m)
                density_estimate=result$dme_un(data)
                true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
                diff_umd=density_estimate-true_density
                sum_e2=0
                for(i in 1:n){
                     sum_e2=sum_e2+(abs(diff_umd[i]))^2
                   }
                RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
         }
         RMISE_estimate=RMISE_SUM/K
         return(RMISE_estimate)}


RMISE_kde2d<- function (K,data){
  RMISE_SUM=0
  for(k in 1:K){
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n^2)
    mass_result=kde2d(x=data[,1],y=data[,2],n = 20, lims = c(range(data[,1]), range(data[,2])))

    # create a new data frame of that 2d density grid
    gr <- data.frame(with(mass_result, expand.grid(x,y)), as.vector(mass_result$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    
    # Apply the model to the original data to estimate density at that point
    smoothed_density<- predict(mod, newdata=data.frame(xgr=data[,1], ygr=data[,2]))
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_kde2d=as.vector(smoothed_density-true_density)
    sum_e2=0
    for(i in 1:n){
      sum_e2=sum_e2+(abs(diff_kde2d[i]))^2
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
  }
  RMISE_estimate=RMISE_SUM/K
  return(RMISE_estimate)}



RMISE_ks<- function (K,data){
  RMISE_SUM=0
  for(k in 1:K){
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    Hpi1 <- Hpi(x=data)
    ks_result= ks::kde(x=data, H=Hpi1,eval.points = data)
    density_estimate=ks_result$estimate
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_ks=density_estimate-true_density
    sum_e2=0
    for(i in 1:n_dim){
      sum_e2=sum_e2+(abs(diff_ks[i]))^2
    }
    RMISE_SUM=RMISE_SUM+sqrt(d*sum_e2)
  }
  RMISE_estimate=RMISE_SUM/K
  return(RMISE_estimate)}

RMISE_kdevine<- function (K,data){
  RMISE_SUM=0
  for(k in 1:K){
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    fit <- kdevine(data)
    density_estimate=dkdevine(data, fit)
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
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
  data=mvrnorm(n=n,mu=c(0,0),Sigma=Sigma)
  rmise3_umd=cbind(rmise3_umd,RMISE_umd(100,data=data))
  rmise3_kde2d=cbind(rmise3_kde2d,RMISE_kde2d(100,data=data))
  rmise3_ks=cbind(rmise3_ks,RMISE_ks(100,data=data))
  rmise3_kdevine=cbind(rmise3_kdevine,RMISE_kdevine(100,data=data))
  
}

df.rmise3_kde2d=data.frame(RMISE=t(rmise3_kde2d),group=rep("kde2d",v0))
df.rmise3_umd=data.frame(RMISE=t(rmise3_umd),group=rep("Berstein",v0))
df.rmise3_ks=data.frame(RMISE=t(rmise3_ks),group=rep("ks",v0))
df.rmise3_kdevine=data.frame(RMISE=t(rmise3_kdevine),group=rep("kdevine",v0))

rmise_gaussian=rbind(df.rmise3_umd,df.rmise3_kde2d,df.rmise3_ks,df.rmise3_kdevine)


#draw boxplot for rmise Gaussian comparision
library(ggplot2)
boxplot_rmise_gaussian=ggplot(rmise_gaussian, aes(x=group,y=RMISE)) + 
                              geom_boxplot()+
                              coord_cartesian(ylim = c(0, max(c(rmise3_umd,rmise3_kde2d,rmise3_ks,rmise3_kdevine))))

png(file = "boxplot_rmise_gaussian.png", bg = "transparent")
boxplot_rmise_gaussian
dev.off()

boxplot(t(rmise3_umd))
boxplot(t(rmise3_kde2d))


########################################################
#The following generates ML1E for Gaussian Distribution#	
###########################################################
ML1E_umd<- function (K,data){
  ML1E_SUM=0
  for(k in 1:K){
    d=((max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2])))/(n^2)
    result=umd_multi(data,crit="CN",m=m)
    density_estimate=result$dme_un(data)
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_umd=density_estimate-true_density
    sum_abe=0
    for(i in 1:n){
      sum_abe=sum_abe+abs(diff_umd[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}


ML1E_kde2d<- function (K,data){
  ML1E_SUM=0
  for(k in 1:K){
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n^2)
    mass_result=kde2d(x=data[,1],y=data[,2],n = 20, lims = c(range(data[,1]), range(data[,2])))
    
    # create a new data frame of that 2d density grid
    gr <- data.frame(with(mass_result, expand.grid(x,y)), as.vector(mass_result$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    
    # Apply the model to the original data to estimate density at that point
    smoothed_density<- predict(mod, newdata=data.frame(xgr=data[,1], ygr=data[,2]))
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_kde2d=as.vector(smoothed_density-true_density)
    sum_abe=0
    for(i in 1:n){
      sum_abe=sum_abe+abs(diff_kde2d[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}


ML1E_ks<- function (K,data){
  ML1E_SUM=0
  for(k in 1:K){
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    Hpi1 <- Hpi(x=data)
    ks_result= ks::kde(x=data, H=Hpi1,eval.points = data)
    density_estimate=ks_result$estimate
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_ks=density_estimate-true_density
    sum_abe=0
    for(i in 1:n_dim){
      sum_abe=sum_abe+abs(diff_ks[i])
    }
    ML1E_SUM=ML1E_SUM+d*sum_abe
  }
  ML1E_estimate=ML1E_SUM/K
  return(ML1E_estimate)}

ML1E_kdevine<- function (K,data){
  ML1E_SUM=0
  for(k in 1:K){
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    fit <- kdevine(data)
    density_estimate=dkdevine(data, fit)
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
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
  data=mvrnorm(n=n,mu=c(0,0),Sigma=Sigma)
  ML1E3_umd=cbind(ML1E3_umd,ML1E_umd(100,data=data))
  ML1E3_kde2d=cbind(ML1E3_kde2d,ML1E_kde2d(100,data=data))
  ML1E3_ks=cbind(ML1E3_ks,ML1E_ks(100,data=data))
  ML1E3_kdevine=cbind(ML1E3_kdevine,ML1E_kdevine(100,data=data))
  
  
}



df.ML1E3_kde2d=data.frame(ML1E=t(ML1E3_kde2d),group=rep("kde2d",v0))
df.ML1E3_umd=data.frame(ML1E=t(ML1E3_umd),group=rep("Berstein",v0))
df.ML1E3_ks=data.frame(ML1E=t(ML1E3_ks),group=rep("ks",v0))
df.ML1E3_kdevine=data.frame(ML1E=t(ML1E3_kdevine),group=rep("kdevine",v0))

ML1E_gaussian=rbind(df.ML1E3_umd,df.ML1E3_kde2d,df.ML1E3_ks,df.ML1E3_kdevine)

#draw boxplot for ML1E Gaussian comparision
library(ggplot2)
boxplot_ML1E_gaussian=ggplot(ML1E_gaussian, aes(x=group,y=ML1E)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML1E3_kde2d,ML1E3_ks,ML1E3_kdevine)))

png(file = "boxplot_ML1E_gaussian.png", bg = "transparent")
boxplot_ML1E_gaussian
dev.off()

boxplot(t(ML1E3_umd))
boxplot(t(ML1E3_kde2d))

########################################################
#The following generates ML_infE for Gaussian Distribution#	
###########################################################
ML_infE_umd<- function (K,data){
  ML_infE_SUM=0
  for(k in 1:K){
    d=((max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2])))/(n^2)
    result=umd_multi(data,crit="CN",m=m)
    density_estimate=result$dme_un(data)
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_umd=density_estimate-true_density
    sum_max_abs=max(abs(diff_umd))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}


ML_infE_kde2d<- function (K,data){
  ML_infE_SUM=0
  for(k in 1:K){
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n^2)
    mass_result=kde2d(x=data[,1],y=data[,2],n = 20, lims = c(range(data[,1]), range(data[,2])))
    
    # create a new data frame of that 2d density grid
    gr <- data.frame(with(mass_result, expand.grid(x,y)), as.vector(mass_result$z))
    names(gr) <- c("xgr", "ygr", "zgr")
    
    # Fit a model
    mod <- loess(zgr~xgr*ygr, data=gr)
    # Apply the model to the original data to estimate density at that point
    smoothed_density<- predict(mod, newdata=data.frame(xgr=data[,1], ygr=data[,2]))
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_kde2d=as.vector(smoothed_density-true_density)
    sum_max_abs=max(abs(diff_kde2d))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}


ML_infE_ks<- function (K,data){
  ML_infE_SUM=0
  for(k in 1:K){
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    Hpi1 <- Hpi(x=data)
    ks_result=ks::kde(x=data, H=Hpi1,eval.points = data)
    density_estimate=ks_result$estimate
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
    diff_ks=density_estimate-true_density
    sum_max_abs=max(abs(diff_ks))
    ML_infE_SUM=ML_infE_SUM+sum_max_abs
  }
  ML_infE_estimate=ML_infE_SUM/K
  return(ML_infE_estimate)}

ML_infE_kdevine<- function (K,data){
  ML_infE_SUM=0
  for(k in 1:K){
    n_dim=nrow(data)
    d=(max(data[,1])-min(data[,1]))*(max(data[,2])-min(data[,2]))/(n_dim^2)
    fit <- kdevine(data)
    density_estimate=dkdevine(data, fit)
    true_density=dmvnorm(data,mean=c(0,0),sigma=Sigma)
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
  data=mvrnorm(n=n,mu=c(0,0),Sigma=Sigma)
  ML_infE3_umd=cbind(ML_infE3_umd,ML_infE_umd(100,data=data))
  ML_infE3_kde2d=cbind(ML_infE3_kde2d,ML_infE_kde2d(100,data=data))
  ML_infE3_ks=cbind(ML_infE3_ks,ML_infE_ks(100,data=data))
  ML_infE3_kdevine=cbind(ML_infE3_kdevine,ML_infE_kdevine(100,data=data))
  
}

df.ML_infE3_kde2d=data.frame(ML_infE=t(ML_infE3_kde2d),group=rep("kde2d",v0))
df.ML_infE3_umd=data.frame(ML_infE=t(ML_infE3_umd),group=rep("Berstein",v0))
df.ML_infE3_ks=data.frame(ML_infE=t(ML_infE3_ks),group=rep("ks",v0))
df.ML_infE3_kdevine=data.frame(ML_infE=t(ML_infE3_kdevine),group=rep("kdevine",v0))

ML_infE_gaussian=rbind(df.ML_infE3_umd,df.ML_infE3_kde2d,df.ML_infE3_ks,df.ML_infE3_kdevine)


#draw boxplot for ML_infE Gaussian comparision
library(ggplot2)
boxplot_ML_infE_gaussian=ggplot(ML_infE_gaussian, aes(x=group,y=ML_infE)) + 
  geom_boxplot()+
  coord_cartesian(ylim = c(0, max(ML_infE3_kde2d,ML_infE3_ks,ML_infE3_kdevine)))

png(file = "boxplot_ML_inf_gaussian.png", bg = "transparent")
boxplot_ML_infE_gaussian
dev.off()

boxplot(t(ML_infE3_umd))
boxplot(t(ML_infE3_kde2d))

##########################################################################
##########################################################################
##########################################################################

# Outlier case
p=rbern(1, prob = 0.95)
Sigma0=matrix(c(1,0,0,1),ncol=2,nrow=2)
set.seed(20240612)
part1=mvrnorm(n=950,mu=c(0,0),Sigma=Sigma0)
part2=mvrnorm(n=50,mu=c(6,6),Sigma=Sigma0)
data_outlier=rbind(part1,part2)
data_outlier=data_outlier[is.finite(data_outlier[,2]) == 1 & is.finite(data_outlier[,1])==1, ]
data_outlier=data.frame(x=data_outlier[,1],y=data_outlier[,2])
plot(x=data_outlier[,1],y=data_outlier[,2])
x=seq(-10,10,0.5)
y=seq(-10,10,0.5)
grid_data<- expand.grid(x,y)
grid_data<-data.frame(x=grid_data[,1],y=grid_data[,2])


n_dim_1=nrow(data_outlier)
d_1=(max(data_outlier[,1])-min(data_outlier[,1]))*(max(data_outlier[,2])-min(data_outlier[,2]))/(n_dim_1^2)

#kdevine
fit_1 <- kdevine(data_outlier)
density_estimate_kedvine_1=dkdevine(grid_data, fit_1)
grid_data$kdevine_estimate=density_estimate_kedvine_1
# true_density=dmvnorm(data_outlier,mean=c(0,0),sigma=Sigma)
# diff_kdevine_1=density_estimate_kedvine-true_density

#ks
Hpi1_1 <- Hpi(x=data_outlier)
ks_result_1=ks::kde(x=data_outlier, H=Hpi1_1,eval.points =grid_data[,1:2])
grid_data$ks_estimate=ks_result_1$estimate
# diff_ks_1=density_estimate_ks-true_density

#kde2d
mass_result_1=kde2d(x=data_outlier[,1],y=data_outlier[,2],n = 20, lims = c(range(data_outlier[,1]), range(data_outlier[,2])))
# create a new data frame of that 2d density grid
gr_1 <- data.frame(with(mass_result_1, expand.grid(x,y)), as.vector(mass_result_1$z))
names(gr_1) <- c("x", "y", "z")
# Fit a model
mod_1 <- loess(z~x*y, data=gr_1)
# Apply the model to the original data to estimate density at that point
smoothed_density_kde2d<- predict(mod_1, newdata=grid_data[,1:2],na.action = na.pass)
smoothed_density_kde2d=smoothed_density_kde2d %>% replace(is.na(.), 0)
# diff_kde2d_1=as.vector(smoothed_density_kde2d-true_density)
grid_data$kde2d_estimate=as.vector(smoothed_density_kde2d)

#umd
result_umd_multi=umd_multi(data_outlier,crit="CN",m=20)
density_estimate_umd_multi=result_umd_multi$dme_un(grid_data[,1:2])
# diff_umd_multi=density_estimate_umd_multi-true_density

grid_data$umd_estimate<-density_estimate_umd_multi

p_umd <-  ggplot() +
  geom_point(data=data_outlier,aes(x,y))+
  geom_contour(data=grid_data, aes(x, y, z = umd_estimate),bins = 20)+  
  scale_x_continuous(expand = c(0, 0),limits = c(-5,5)) +
  scale_y_continuous(expand = c(0, 0),limits = c(-5,5)) +# Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  ggtitle('Berstein')+
  theme_bw() +                                                      # dark-on-light theme
  theme(panel.background = element_blank(),                    
        axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        plot.title = element_text(face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))

p_umd

p_kde2d<- ggplot() +
  geom_point(data=data_outlier,aes(x,y))+
  geom_contour(data=grid_data,aes(x, y, z =kde2d_estimate),bins = 20)+               # Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  ggtitle('kde2d')+
  scale_x_continuous(expand = c(0, 0),limits = c(-5,5)) +
  scale_y_continuous(expand = c(0, 0),limits = c(-5,5)) +# Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  theme_bw() +                                                      # dark-on-light theme
  theme(panel.background = element_blank(),                    
        axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        plot.title = element_text(face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))

p_kde2d

p_ks<- ggplot() +
  geom_point(data=data_outlier,aes(x,y))+
  geom_contour(data=grid_data,aes(x, y, z=ks_estimate),bins = 20)+               # Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  ggtitle('ks')+
  scale_x_continuous(expand = c(0, 0),limits = c(-5,5)) +
  scale_y_continuous(expand = c(0, 0),limits = c(-5,5)) +# Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  theme_bw() +                                                      # dark-on-light theme
  theme(panel.background = element_blank(),                    
        axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        plot.title = element_text(face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))

p_ks

p_kdevine<- ggplot() +
  geom_point(data=data_outlier,aes(x,y))+
  geom_contour(data=grid_data,aes(x, y, z=kdevine_estimate),bins = 20)+               # Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  ggtitle('kdevine')+
  scale_x_continuous(expand = c(0, 0),limits = c(-5,5)) +
  scale_y_continuous(expand = c(0, 0),limits = c(-5,5)) +# Create contour plot with 15 bins          # create contour plot with 30 bins and set color of level
  theme_bw() +                                                      # dark-on-light theme
  theme(panel.background = element_blank(),                    
        axis.line.x = element_line(),                               # these two are for the axis line
        axis.line.y = element_line(),
        axis.text.x = element_text(colour = "black"),               # there two are for texts in axes
        axis.text.y = element_text(colour = "black"),
        axis.ticks.x = element_line(),                              # these two are for ticks in axes
        axis.ticks.y = element_line(),
        axis.title.x = element_text(colour = "black", face = 'bold', vjust = -1),                              
        axis.title.y = element_text(colour = "black", face = 'bold'),
        plot.title = element_text(face = 'bold'),
        legend.title = element_text(colour = "black", face = 'bold'),
        legend.text = element_text(colour = "black"))

p_kdevine


library(ggExtra)
marginal_plot_kde2d=ggMarginal(p_kde2d)
marginal_plot_kde2d
marginal_plot_umd=ggMarginal(p_umd)
marginal_plot_umd
marginal_plot_ks=ggMarginal(p_ks)
marginal_plot_ks
marginal_plot_kdevine=ggMarginal(p_kdevine)
marginal_plot_kdevine

marginalx_grid_data=grid_data%>% group_by(x) %>%
  summarise(marginalx_kde2d= sum(kde2d_estimate),
            marginalx_ks= sum(ks_estimate),
            marginalx_umd= sum(umd_estimate),
            marginalx_kdevine= sum(kdevine_estimate))
#melt data frame into long format
df_marginal_x<- melt(marginalx_grid_data , id.vars = 'x', variable.name = 'Method')

#create line plot for each column in data frame
png("X_marginal_density.png", width=3000, height=2000, res=800, units="px")
ggplot(df_marginal_x, aes(x, value)) +
  geom_line(aes(colour = Method))+
  theme_classic()+
  ylab("Marginal density of X ")+
  scale_color_manual(name='Method', labels=c('kde2d', 'ks', 'Bernstein', 'kdevine'), values=c('black', 'purple', 'red','blue'))
dev.off()

save.image("20240610_gaussian_outlier_seed.RData")
