
###############################################################################
#            Multivariate Density Estimation Using Bernstein Polynomials          #
#																	      	  #
#          Copyright(2021): Dan Han, University of Louisville, KY, USA
#																			  #	
#                         Version 2 - June 28, 2021					      #	
###############################################################################

#Packages needed
require(quadprog)
require(bivariate)
require(magic)

#Function to determine the optimal number of weights
#using the CN criterion
mOpt.CN <- function(data, D){
  
  #set starting point for m
  n <- nrow(data)
  m <- floor(  n^(1/3) ) - 1
  
  sort.data1 <- sort(data[,1])
  sort.data2 <- sort(data[,2])
  
  #Make the transformed data set
  tdata1 <- (data[,1]-sort.data1[1])/(sort.data1[n]-sort.data1[1])	
  tdata2 <- (data[,2]-sort.data2[1])/(sort.data2[n]-sort.data2[1])	
  #set start value for while loop
  logratio <- 1
  
  while( logratio < sqrt(n) ){
    m <- m+1
    
    
    #construct A matrix using n value 
    G <- NULL
    for(i in 1:n){
      F <-NULL
      for(k1 in 0:m){
        for(k2 in 0:m){
        F<-rbind(F,pbeta(tdata1[i], shape1=k1+1, shape2=m-k1+1)*pbeta(tdata2[i], shape1=k2+1, shape2=m-k2+1))
        }
      }
      G <- cbind(G, F)
      
    }
    G=t(G)
    Amat <- t(G) %*% D %*% G		
    
    #take spectral decomposition to find eigenvalues
    spec <- eigen(Amat, symmetric=TRUE)
    d <- spec$values
    min.eigenValue <- max( min(d), 0 )
    max.eigenValue <- max(d)
    
    logratio <- log10(max.eigenValue) - log10(min.eigenValue)		
  }	
  m-1 #return number of weights
}



#Function to determine the index of the maximum weight
maxWeight <- function(m, Fn){
  #find max of weights
  maxplace <- which.max( Fn( ((1:m)/m)) - 
                           Fn( ((0:(m-1))/m) )  )
}

########################### wrong###########################
#Function to generate the constraint matrix 
constraintMat <- function( m, maxplace_1,maxplace_2){
  Q1=kronecker(diag(m+1),t(rep(1,m+1)))
  Q2=do.call("cbind", rep(list(diag(m+1)), m+1))
  C1=suppressWarnings(rbind(matrix( rep( c(-1,1, rep(0,m)) , maxplace_1), maxplace_1, m+1, byrow=TRUE), 
            matrix( rep( c( rep (0,maxplace_1),1,-1,rep(0, m-maxplace_1)),m-maxplace_1),
                    m-maxplace_1,m+1,byrow=TRUE)))
  C2=suppressWarnings(rbind(matrix( rep( c(-1,1, rep(0,m)) , maxplace_2), maxplace_2, m+1, byrow=TRUE), 
                            matrix( rep( c( rep (0,maxplace_2),1,-1,rep(0, m-maxplace_2)),m-maxplace_2),
                                    m-maxplace_2,m+1,byrow=TRUE)))
  # the following is R matrix in the paper
  R <-    rbind( rep(1,(m+1)^2), diag(rep(1,(m+1)^2)), 
          C1 %*% Q1,
          C2 %*% Q2)
  #the constraint matrix in solve.QP function, it is t(A)w >= w_0, t(A) is R here, thus the constraint matrix A in solve.QP is t(R)=Rmat
  Rmat <- t(R)
}


#Function to solve for the weight vector	
solveWeights <- function(m, Fn1, Fn2, Amat, bvec){	
  #find the location of the maximum weight
  max.place1 <- maxWeight(m, Fn1)
  max.place2<-maxWeight(m, Fn2)
  #make the constraint matrix
  Rmat <- constraintMat(m, max.place1,max.place2)	
  #make rvec vector of constraints
  rvec=c(1,rep(0,(m+1)^2+2*m))
  #Find Dmat matrix for solve.QP function
  Dmat=2*Amat
  #Find dvec matrix for solve.QP function
  dvec=t(bvec)
  #find weights using solve.QP function
  w.hat = solve.QP(Dmat,dvec,Rmat,rvec,meq=1)$solution
  #function to find max of an element and 0
  max0 <- function(x){max(x,0)}
  #make sure no weights are < 0
  w.hat <- sapply( w.hat, max0)
  #make sure the weights sum to 1
  wsum <- sum(w.hat)
  w.hat <- w.hat / wsum
}	



#The main function which the user will call
umd_multi <- function(data,crit="CN", m=NA, warning=TRUE){

  #delta definition
  n <- nrow(data)
  delta1 <- sd(data[,1])/ sqrt(n)
  delta2 <- sd(data[,2])/ sqrt(n)
  sort.data1 <- sort(data[,1])
  sort.data2 <- sort(data[,2])
  
  #Make the transformed data set
  tdata1 <- (data[,1]-sort.data1[1])/(sort.data1[n]-sort.data1[1])	
  tdata2 <- (data[,2]-sort.data2[1])/(sort.data2[n]-sort.data2[1])	
  
  #Construct the D matrix and make vector of empirical cdf values
  Fn<-ebvcdf(x=tdata1,y=tdata2)  #Fn is empirical cdf 
  ep <- 3/( 8*n )	
  ecdf.vec <- Fn(tdata1,tdata2)	
  g <- ecdf.vec
  
  #construct D
  D <- diag( n / ( ( ecdf.vec + ep)*(1 + ep - ecdf.vec ) ) )
  
  
  #Find the optimal number of weights, depends on the given crit
  #then return the desired functions
  if( crit == "CN"){
    mOpt <- mOpt.CN(data, D)
    
    if( is.na(m) == FALSE){
      
      if( mOpt < m && warning==TRUE){
        cat("WARNING: given number of weights is larger than optimal number,\n"
            , "\t \t optimal number of weights =", mOpt, "\n")
      }
      if( mOpt > m && warning==TRUE){
        cat("WARNING: given number of weights is less than optimal number,\n"
            ,"\t \t optimal number of weights =", mOpt, "\n")
      }
      
      #construct A matrix using n value 
      G <- NULL
      for(i in 1:n){
        F <-NULL
        for(k1 in 0:m){
          for(k2 in 0:m){
            F<-rbind(F,pbeta(tdata1[i], shape1=k1+1, shape2=m-k1+1)*pbeta(tdata2[i], shape1=k2+1, shape2=m-k2+1))
          }
        }
        G <- cbind(G, F)
        
      }
      G=t(G)
      Amat <- t(G) %*% D %*% G				
      bvec <- 2*t(g)%*% D %*% G
      
      #Make sure there are no zero eigen values	
      spec <- eigen(Amat, symmetric=TRUE)
      d <- spec$values
      Q <- spec$vectors
      
      #find which values are < 10e-6, and set to smallest
      #eigen value >= 10e-6
      if( min(d) < 10e-6){
        tooSmall <- which( d < 10e-6)	
        d[tooSmall] <- d[ min(tooSmall) - 1]
        #Recreate pos. def. Amat matrix 
        Amat <- Q %*% diag(d) %*% t(Q)				
      }	
      
    } #if no m provided then use the optimal m
    else{ m <- mOpt 	
    
    #make G, Amat, and bvec
    #construct A matrix using n value 
    G <- NULL
    for(i in 1:n){
      F <-NULL
      for(k1 in 0:m){
        for(k2 in 0:m){
          F<-rbind(F,pbeta(tdata1[i], shape1=k1+1, shape2=m-k1+1)*pbeta(tdata2[i], shape1=k2+1, shape2=m-k2+1))
        }
      }
      G <- cbind(G, F)
      
    }
    G=t(G)
    Amat <- t(G) %*% D %*% G				
    bvec <- 2*t(g)%*% D %*% G
    
    #Make sure there are no zero eigen values	
    spec <- eigen(Amat, symmetric=TRUE)
    d <- spec$values
    Q <- spec$vectors
    
    #find which values are < 10e-6, and set to smallest
    #eigen value >= 10e-6
    if( min(d) < 10e-6){
      tooSmall <- which( d < 10e-6)	
      d[tooSmall] <- d[ min(tooSmall) - 1]
      #Recreate pos. def. Amat matrix 
      Amat <- Q %*% diag(d) %*% t(Q)
    }
    }
    
    #Solve for the weights
    Fn1=ecdf(tdata1)
    Fn1_val=Fn1(tdata1)
    Fn2=ecdf(tdata2)
    Fn2_val=Fn2(tdata2)
    
    weights <- solveWeights(m,Fn1,Fn2,Amat,bvec)

    #use the weights to create the 4 distribution functions
    #dme_un stands for multivariate density estimation subject to Marginal Unimodality Constraints
    dme_un=function(x){
      
      mix.pdf=function(x){
      B_m <- NULL
      for(k1 in 0:m){
        for(k2 in 0:m){
            B_m<-rbind(B_m,dbeta((x[1]-sort.data1[1])/(sort.data1[n]-sort.data1[1]), shape1=k1+1, shape2=m-k1+1)*dbeta((x[2]-sort.data2[1])/(sort.data2[n]-sort.data2[1]), shape1=k2+1, shape2=m-k2+1))
        }
      }
        density_umd=weights%*%B_m/((sort.data1[n]-sort.data1[1])*(sort.data2[n]-sort.data2[1]))
          }
  
      apply(x,1, mix.pdf)
      }
      
    
    
    return(list(weights=weights,m.hat=m,dme_un=dme_un))
  }}

