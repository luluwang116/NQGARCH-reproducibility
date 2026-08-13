
Wmatrix0 <- read.csv( "simulateweight.csv",header = TRUE)

Wmatrix <- as.matrix(Wmatrix0[, -1])

storage.mode(Wmatrix) <- "double"


make_sparse_W <- function(W, K = 3){
  
  W_sparse <- as.matrix(W)
  
  diag(W_sparse) <- 0
  
  N <- nrow(W_sparse)
  
  for(i in 1:N){
    
    candidates <- setdiff(1:N, i)
    
    ord <- candidates[
      order(
        W_sparse[i, candidates],
        decreasing = TRUE
      )
    ]
    
    keep <- ord[
      1:min(K, length(ord))
    ]
    
    drop <- setdiff(
      1:N,
      c(i, keep)
    )
    
    W_sparse[i, drop] <- 0
    W_sparse[i, i] <- 0
  }
  
  rs <- rowSums(W_sparse)
  
  W_sparse <- W_sparse / rs
  
  return(W_sparse)
}


W <- make_sparse_W( Wmatrix, K = 3)

garchcoef <- read.csv(
  "simulation_coef_network.csv",
  header = TRUE
)

garchcoef <- as.matrix(garchcoef)

storage.mode(garchcoef) <- "double"


alpha0 <- garchcoef[, 1]

alpha1 <- garchcoef[, 2]

beta <- garchcoef[,3]

rho <- garchcoef[,4]

h0<-alpha0/(1-alpha1-beta-rho)

###########################################################
set.seed(17234567)

n <- 30
T <- 1050
S <- 100

xit1 <- array(0, dim = c(n, T, S))
eta1 <- array(0, dim = c(n, T, S))
h1   <- array(0, dim = c(n, T, S))

############################################################

eta1[] <- rnorm(
  length(eta1),0,1)

############################################################
# Generate Network-GARCH process
############################################################

for(s in 1:S){
  
  h1[, 1, s] <-
    alpha0 +
    alpha1 * h0 +
    beta * h0 +
    rho * as.numeric(
      W %*% h0
    )
  
  xit1[, 1, s] <-
    sqrt(h1[, 1, s]) *
    eta1[, 1, s]
  
  
  for(t in 2:T){
    
    x2_lag <- xit1[, t - 1, s]^2
    
    h1[, t, s] <-
      alpha0 +
      alpha1 * x2_lag +
      beta * h1[, t - 1, s] +
      rho * as.numeric(
        W %*% x2_lag
      )
    
    xit1[, t, s] <-
      sqrt(h1[, t, s]) *
      eta1[, t, s]
  }
}


############################################################
# Remove first 50 observations
############################################################

xit1 <- xit1[, 51:T, , drop = FALSE]
h1   <- h1[, 51:T, , drop = FALSE]
eta1 <- eta1[, 51:T, , drop = FALSE]


dim(xit1)
dim(h1)
dim(eta1)
########################################
set.seed(29234567)

n <- 30
T <- 1050
S <- 100

xit2 <- array(0, dim = c(n, T, S))
eta2 <- array(0, dim = c(n, T, S))
h2   <- array(0, dim = c(n, T, S))


eta2[] <- rnorm(
  length(eta2),0,1)

############################################################
# Generate Network-GARCH process
############################################################

for(s in 1:S){
  
  h2[, 1, s] <-
    alpha0 +
    alpha1 * h0 +
    beta * h0 +
    rho * as.numeric(
      W %*% h0
    )
  
  xit2[, 1, s] <-
    sqrt(h2[, 1, s]) *
    eta2[, 1, s]
  
  
  for(t in 2:T){
    
    x2_lag <- xit2[, t - 1, s]^2
    
    h2[, t, s] <-
      alpha0 +
      alpha1 * x2_lag +
      beta * h2[, t - 1, s] +
      rho * as.numeric(
        W %*% x2_lag
      )
    
    xit2[, t, s] <-
      sqrt(h2[, t, s]) *
      eta2[, t, s]
  }
}


############################################################
# Remove first 50 observations
############################################################

xit2 <- xit2[, 51:T, , drop = FALSE]
h2   <- h2[, 51:T, , drop = FALSE]
eta2 <- eta2[, 51:T, , drop = FALSE]


dim(xit2)
dim(h2)
dim(eta2)
##############################################
library(abind)
x<-abind(xit1,xit2,along = 3)
h<-abind(h1,h2,along=3)
dim(x)

truequantile01<-qnorm(0.1)*sqrt(h)

truequantile005<-qnorm(0.05)*sqrt(h)
###################NQGARCH-FORECAST######################
library(quantreg)
library(rugarch)
######################
h_update<-function(N,T,W,y,theta,intc)  #iterative formula of h_t for GARCH(1,1)
{
  #N is spatial sample size;T is time series sample size;
  #W is weighted matrix;y is spatio-temporal data from GARCH(1,1); 
  #theta is the parameter; (p+q+1+1)*1 vector 
  #intc is the initial value (N*T)
  alpha0=theta[,1];alpha=theta[,2];beta=theta[,3];rho=theta[,4]
  sigma=matrix(NA,N,T)
  for (i in 1:N) {
    sigma[i,1]=alpha0[i]+alpha[i]*intc[i]+beta[i]*intc[i]+rho[i]*(W[i,]%*%as.matrix(intc,N,1))
    #initial values for y_0^2 and h_0, both are intc for each spatial index i
  }
  
  for (t in 2:T)
  {
    for (i in 1:N) {
      sigma[i,t]=alpha0[i]+alpha[i]*y[i,t-1]^2+beta[i]*sigma[i,t-1]+rho[i]*(W[i,]%*%(y[,t-1]^2))
    }
    
  } 
  return(sigma) #h_t
}
h_update1<-function(T,Wi,xi,X,thetai,intci,intc)  #iterative formula of h_t for GARCH(1,1)
{
  #N is spatial sample size;T is time series sample size;
  #W is weighted matrix;y is spatio-temporal data from GARCH(1,1); 
  #theta is the parameter; (p+q+1+1)*1 vector 
  #intc is the initial value (N*T)
  alpha0=thetai[1];alpha=thetai[2];beta=thetai[3];rho=thetai[4]
  sigma1=matrix(NA,1,T)
  sigma1[1,1]=alpha0+alpha*intci+beta*intci+rho*(Wi%*%as.matrix(intc,N,1))
  #initial values for y_0^2 and h_0, both are intc for each spatial index i
  
  for (t in 2:T)
  {
    sigma1[1,t]=alpha0+alpha*xi[t-1]^2+beta*sigma1[1,t-1]+rho*(Wi%*%(X[,t-1]^2))
    
  } 
  return(sigma1) #h_t
}
######################
derivativeG<-function(N,T,W,x,theta,intc)  #score vector and Hessian for Gaussian GARCH(1,1)
{
  #N is sample size; orders p,q=1; x is N*T matrix
  #theta is the QMLE parameter; (p+q+1+1)*1 vector  
  htilde=h_update(N,T,W,x,theta,intc) #N*T vector
  X1=cbind(intc,x[,1:(T-1)]^2) #initial value for x_0^2 is intc
  X2=cbind(intc,htilde[,1:(T-1)]) #initial value for htilde_0 is intc
  X3=W%*%(cbind(intc,x[,1:(T-1)]^2))#initial value for spatial lag
  X1C=as.matrix(X1)
  X2C=as.matrix(X2)
  X3C=as.matrix(X3)
  beta=theta[,3]
  
  P1h=array(0,dim=c(N,T,4)) #first derivative of h_t with respect to theta given initial values of P1h_0 being zero
  P2h=array(0,dim=c(N,T,4,4)) #second derivative of h_t with respect to theta given initial values of P2h_0 being zero 
  P1ell=array(0,dim=c(N,T,4)) #first derivative of ell_t with respect to theta (score vector)
  P2ell=array(0,dim=c(N,T,4,4)) #second derivative of ell_t with respect to theta (Hessian)
  for (i in 1:N) {
    P1h[i,1,]=c(1,X1[i,1],X2[i,1],X3[i,1])
  }
  for (i in 1:N) {
    P1ell[i,1,]=(1-x[i,1]^2/htilde[i,1])/htilde[i,1]*P1h[i,1,]
    P2ell[i,1,,]=(1-x[i,1]^2/htilde[i,1])/htilde[i,1]*P2h[i,1,,]+(2*x[i,1]^2/htilde[i,1]-1)/(htilde[i,1]^2)*tcrossprod(P1h[i,1,])
  }
  
  
  for (t in 2:T)
  {   
    
    for (i in 1:N) {
      P1h[i,t,]=c(1,X1[i,t],X2[i,t],X3[i,t])+beta[i]*P1h[i,t-1,]
      P2h[i,t,,]=cbind(0,0,P1h[i,t-1,],0)+rbind(0,0,P1h[i,t-1,],0)+beta[i]*P2h[i,t-1,,]    
      
      P1ell[i,t,]=(1-x[i,t]^2/htilde[i,t])/htilde[i,t]*P1h[i,t,]
      P2ell[i,t,,]=(1-x[i,t]^2/htilde[i,t])/htilde[i,t]*P2h[i,t,,]+(2*x[i,t]^2/htilde[i,t]-1)/(htilde[i,t]^2)*tcrossprod(P1h[i,t,])
    }
  }
  
  return(list(P1h=P1h,P2h=P2h,P1ell=P1ell,P2ell=P2ell))
} 
###################################
QMLE_NR<-function(N,T,W,MM,p,q,tolerance,x,theta_int,intc)  #QMLE for GARCH(1,1)
{
  #N is sample size; orders p,q=1
  #MM is maximum iteration
  #tolerance is the convergence tolerance
  #x is t.s. data; N*T matrix   
  #theta_int is initial values for parameter; (p+q+1+1)*1 vector
  
  theta=array(0,dim=c(N,p+q+1+1,MM))
  theta_up=matrix(0,N,4)
  Indicator=rep(0,N)
  laststep=rep(1,N) 
  
  pb <- txtProgressBar(min = 0, max = N,style = 3, width = 50, char = "=")
  for (i in 1:N) {
    s=1 
    theta[,,s]<-theta_int
    repeat
    {
      der=derivativeG(N,T,W,x,theta[,,s],intc)
      P1ells=apply(der$P1ell,c(1,3),sum) #summation of P1ell respect to time t 
      P2ells=apply(der$P2ell,c(1,3,4),sum) #summation of P2ell respect to time t 
      
      
      if (rcond(P2ells[i,,])<10^(-15))
      {
        #print("P2ells[i,,] is computationally singular")
        Indicator[i]=1
        break
      }
      
      theta[i,,s+1]=theta[i,,s]-solve(P2ells[i,,], tol = 1e-16)%*%P1ells[i,]
      
      if (sum(theta[i,,s+1]>0)<4 | sum(theta[i,2:4,s+1]>1)>=1)
      {
        #print("The parameter has values being negative or greater than one")
        Indicator[i]=1
        break    
      } #if the parameter has values being negative or greater than one, stop the repetition.
      else
      {
        laststep[i]=s+1
        maxabs=max(abs(theta[i,,s+1]-theta[i,,s]))
        s=s+1
        #print(S)
        if (s>=MM & maxabs>tolerance)
        {
          Indicator[i]=2
          break
        } #if the iteration approaches the maximum iteration limitation but not converges, stop the repetition.
        if (maxabs<=tolerance) 
        { 
          Indicator[i]=0
          break 
        } #if the iteration converges, stop the repetition.                 
      }    
    }
    
    theta_up[i,]<-theta[i,,laststep[i]] 
    setTxtProgressBar(pb, i)
  }
  return(list(Indicator=Indicator,theta_up=theta_up,laststep=laststep)) 
} 

#####################################
psi<-function(x,tau)
{
  tau-1*(x<0)
}
rk<-function(T,k,ehat0,ehat,tau,weight)
{ 
  #tau is a scalar; k is the lag
  #ehat0 are residuals from estimation; ehat are residuals or bootstrap residuals   
  #weight: (length(ehat)-k)*1 vector of equal weight in diagnosis or the random weight in bootstrap 
  var_a=var(abs(ehat0))
  rk=sum(weight*psi(ehat[(k+1):T],tau)*abs(ehat[1:(T-k)]))/(T)/sqrt((tau-tau^2)*var_a)
  
  return(rk)
}
#########################################
CICQE<-function(Y)
{
  M<-100
  NN<-900
  a<-1
  b<-1000
  N<-30
  tau1<-0.1
  tau2<-0.05
  p<-1
  q<-1
  B<-1000
  
  thetaH1<-thetaH2<-array(0,dim=c(N,p+q+1+1,M)) #hybrid and QMLE
  QuantileH1<-QuantileH2<-matrix(0,N,M) #one-step ahead CQE and its CI 
  #one-step ahead CQE in bootstrap
  QHY1<-QHY2<-matrix(0,N,M)
  
  NN=b-M  #start point of prediction in original data
  n0=max(p,q)
  Ind<-matrix(0,N,M)
  
  for (l in 1:M)
  {
    Xlog=Y[,a:(NN+l-1)]
    T=dim(Xlog)[2]
    X<-Xlog
    Theta_int<-matrix(0,N,p+q+1+1)
    ##############################Data generation############################ 
    for (i in 1:N) {
      #sample size
      ##############################Estimation#################################
      ########################GARCH(1,1) hybrid estimation#####################
      ####step I. QMLE for GARCH(1,1) to obtain h_t
      spatial_lag<-t(W%*%(X^2))[-T,]
      spec1_1=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(p, q),external.regressors=as.matrix(spatial_lag[,i])),                 
                         mean.model = list(armaOrder = c(0, 0),include.mean = FALSE),                 
                         distribution.model = "norm") 
      XG<-X[,-1]
      garch1_1=ugarchfit(spec = spec1_1, data = XG[i,], solver = 'hybrid', fit.control = list(stationarity = 1))
      Theta_int[i,]=as.vector(coef(garch1_1)) #initial value for QMLE
      intc=rowMeans(X[,1:5]^2) #initial value of X_t^2 and h_t for t<=0 
      #initial value can also be mean(X^2) or X[1]^2 or Theta_int[1]/(1-Theta_int[2]-Theta_int[3]) 
      #the influence of intial value is small
    }
    QMLE=QMLE_NR(N,T,W,MM=20,p,q,tolerance=10^(-6),X,Theta_int,intc)    
    QMLEInd=QMLE$Indicator
    thetatilde<-matrix(0,N,4)
    
    XX=array(0,dim=c(N,T,3))
    YY<-matrix(0,N,T)

    for (i in 1:N) {
      
      if(QMLEInd[i]!=0)
      {          
        Ind[i,l]=1#;print(Ind)
        thetatilde[i,]=Theta_int[i,]  #(p+q+1)*1 vector
      }else
      {
        Ind[i,l]=0#;print(Ind)
        thetatilde[i,]=QMLE$theta_up[i,]  #(p+q+1)*1 vector
      }    
      httilde=h_update(N,T,W,X,thetatilde,intc)  #length(httilde)=N  
      
      ####step II. WCQE at specific tau
      #Create Responose and Regressor 
      X1=cbind(intc,X[,1:(T-1)]^2)
      X2=cbind(intc,httilde[,1:(T-1)])
      X3=W%*%X1
      #Regressor without intercept in step II of GARCH(1,1) model; n*(p+q)  
      
      XX[i,,]=cbind(X1[i,],X2[i,],X3[i,])
      YY[i,]=X[i,]^2*sign(X[i,]) #Responose
      fit1=rq(YY[i,] ~ XX[i,,], tau=tau1, weights=1/httilde[i,])
      fit2=rq(YY[i,] ~ XX[i,,], tau=tau2, weights=1/httilde[i,])
      
      thetaH1[i,,l]=as.vector(coef(fit1)) 
      thetaH2[i,,l]=as.vector(coef(fit2)) 
      #############################One-step ahead CQE############################
      QHY1[i,l]=as.numeric(crossprod(thetaH1[i,,l],c(1,X[i,T]^2,httilde[i,T],W[i,]%*%(X[,T]^2))))
      QuantileH1[i,l]=sqrt(abs(QHY1[i,l]))*sign(QHY1[i,l])
      QHY2[i,l]=as.numeric(crossprod(thetaH2[i,,l],c(1,X[i,T]^2,httilde[i,T],W[i,]%*%(X[,T]^2))))
      QuantileH2[i,l]=sqrt(abs(QHY2[i,l]))*sign(QHY2[i,l])
    }
  }
  return(list(QuantileH1=QuantileH1,QuantileH2=QuantileH2))
}
###############################################################
################################
library(parallel)
ncores<-100
index_chunks<-split(1:200,cut(1:200,ncores,labels = FALSE))
cl<-makeCluster(ncores)
clusterEvalQ(cl,{library(quantreg)
  library(rugarch)})
clusterExport(cl,c("CICQE", "h_update", "h_update1", "derivativeG", "QMLE_NR","psi", "rk", "x", "W"))
results01<-parLapply(cl,index_chunks,function(idx){lapply(idx,function(i)CICQE(x[,,i]))})
stopCluster(cl)
##############################################
results0101<-array(0,dim = c(30,100,200))
for (i in 1:100) {
  for (j in 1:2) {
    results0101[,,2*(i-1)+j]<-results01[[i]][[j]]$QuantileH1
  }
}
##############################################
results01005<-array(0,dim = c(30,100,200))
for (i in 1:100) {
  for (j in 1:2) {
    results01005[,,2*(i-1)+j]<-results01[[i]][[j]]$QuantileH2
  }
}
################################################
load("raw_results/Table1/results0101.RData")
load("raw_results/Table1/results01005.RData")
mae1<-mae2<-rep(0,30)
for (i in 1:30) {
  mae1[i]<-mean(abs(results01005[i,,]-truequantile005[i,901:1000,]))
}
for (i in 1:30) {
  mae2[i]<-mean(abs(results0101[i,,]-truequantile01[i,901:1000,]))
}
mae1
mae2
#########################################################
############################QGARCH###################################

library(quantreg)
library(rugarch)
##############################
h_update<-function(N,T,y,theta,intc)  #iterative formula of h_t for GARCH(1,1)
{
  #N is spatial sample size;T is time series sample size;
  #W is weighted matrix;y is spatio-temporal data from GARCH(1,1); 
  #theta is the parameter; (p+q+1+1)*1 vector 
  #intc is the initial value (N*T)
  alpha0=theta[,1];alpha=theta[,2];beta=theta[,3]
  sigma=matrix(NA,N,T)
  for (i in 1:N) {
    sigma[i,1]=alpha0[i]+alpha[i]*intc[i]+beta[i]*intc[i]
    #initial values for y_0^2 and h_0, both are intc for each spatial index i
  }
  
  for (t in 2:T)
  {
    for (i in 1:N) {
      sigma[i,t]=alpha0[i]+alpha[i]*y[i,t-1]^2+beta[i]*sigma[i,t-1]
    }
    
  } 
  return(sigma) #h_t
}
h_update1<-function(T,xi,X,thetai,intci,intc)  #iterative formula of h_t for GARCH(1,1)
{
  #N is spatial sample size;T is time series sample size;
  #W is weighted matrix;y is spatio-temporal data from GARCH(1,1); 
  #theta is the parameter; (p+q+1+1)*1 vector 
  #intc is the initial value (N*T)
  alpha0=thetai[1];alpha=thetai[2];beta=thetai[3]
  sigma1=matrix(NA,1,T)
  sigma1[1,1]=alpha0+alpha*intci+beta*intci
  #initial values for y_0^2 and h_0, both are intc for each spatial index i
  
  for (t in 2:T)
  {
    sigma1[1,t]=alpha0+alpha*xi[t-1]^2+beta*sigma1[1,t-1]
    
  } 
  return(sigma1) #h_t
}
######################
derivativeG<-function(N,T,x,theta,intc)  #score vector and Hessian for Gaussian GARCH(1,1)
{
  #N is sample size; orders p,q=1; x is N*T matrix
  #theta is the QMLE parameter; (p+q+1+1)*1 vector  
  htilde=h_update(N,T,x,theta,intc) #N*T vector
  X1=cbind(intc,x[,1:(T-1)]^2) #initial value for x_0^2 is intc
  X2=cbind(intc,htilde[,1:(T-1)]) #initial value for htilde_0 is intc
  X1C=as.matrix(X1)
  X2C=as.matrix(X2)
  beta=theta[,3]
  
  P1h=array(0,dim=c(N,T,3)) #first derivative of h_t with respect to theta given initial values of P1h_0 being zero
  P2h=array(0,dim=c(N,T,3,3)) #second derivative of h_t with respect to theta given initial values of P2h_0 being zero 
  P1ell=array(0,dim=c(N,T,3)) #first derivative of ell_t with respect to theta (score vector)
  P2ell=array(0,dim=c(N,T,3,3)) #second derivative of ell_t with respect to theta (Hessian)
  for (i in 1:N) {
    P1h[i,1,]=c(1,X1[i,1],X2[i,1])
  }
  for (i in 1:N) {
    P1ell[i,1,]=(1-x[i,1]^2/htilde[i,1])/htilde[i,1]*P1h[i,1,]
    P2ell[i,1,,]=(1-x[i,1]^2/htilde[i,1])/htilde[i,1]*P2h[i,1,,]+(2*x[i,1]^2/htilde[i,1]-1)/(htilde[i,1]^2)*tcrossprod(P1h[i,1,])
  }
  
  
  for (t in 2:T)
  {   
    
    for (i in 1:N) {
      P1h[i,t,]=c(1,X1[i,t],X2[i,t])+beta[i]*P1h[i,t-1,]
      P2h[i,t,,]=cbind(0,0,P1h[i,t-1,])+rbind(0,0,P1h[i,t-1,])+beta[i]*P2h[i,t-1,,]    
      
      P1ell[i,t,]=(1-x[i,t]^2/htilde[i,t])/htilde[i,t]*P1h[i,t,]
      P2ell[i,t,,]=(1-x[i,t]^2/htilde[i,t])/htilde[i,t]*P2h[i,t,,]+(2*x[i,t]^2/htilde[i,t]-1)/(htilde[i,t]^2)*tcrossprod(P1h[i,t,])
    }
  }
  
  return(list(P1h=P1h,P2h=P2h,P1ell=P1ell,P2ell=P2ell))
} 
###################################
QMLE_NR<-function(N,T,MM,p,q,tolerance,x,theta_int,intc)  #QMLE for GARCH(1,1)
{
  #N is sample size; orders p,q=1
  #MM is maximum iteration
  #tolerance is the convergence tolerance
  #x is t.s. data; N*T matrix   
  #theta_int is initial values for parameter; (p+q+1+1)*1 vector
  
  theta=array(0,dim=c(N,p+q+1,MM))
  theta_up=matrix(0,N,3)
  Indicator=rep(0,N)
  laststep=rep(1,N) 
  
  pb <- txtProgressBar(min = 0, max = N,style = 3, width = 50, char = "=")
  for (i in 1:N) {
    s=1 
    theta[,,s]<-theta_int
    repeat
    {
      der=derivativeG(N,T,x,theta[,,s],intc)
      P1ells=apply(der$P1ell,c(1,3),sum) #summation of P1ell respect to time t 
      P2ells=apply(der$P2ell,c(1,3,4),sum) #summation of P2ell respect to time t 
      
      
      if (rcond(P2ells[i,,])<10^(-15))
      {
        #print("P2ells[i,,] is computationally singular")
        Indicator[i]=1
        break
      }
      
      theta[i,,s+1]=theta[i,,s]-solve(P2ells[i,,], tol = 1e-16)%*%P1ells[i,]
      
      if (sum(theta[i,,s+1]>0)<3 | sum(theta[i,2:3,s+1]>1)>=1)
      {
        #print("The parameter has values being negative or greater than one")
        Indicator[i]=1
        break    
      } #if the parameter has values being negative or greater than one, stop the repetition.
      else
      {
        laststep[i]=s+1
        maxabs=max(abs(theta[i,,s+1]-theta[i,,s]))
        s=s+1
        #print(S)
        if (s>=MM & maxabs>tolerance)
        {
          Indicator[i]=2
          break
        } #if the iteration approaches the maximum iteration limitation but not converges, stop the repetition.
        if (maxabs<=tolerance) 
        { 
          Indicator[i]=0
          break 
        } #if the iteration converges, stop the repetition.                 
      }    
    }
    
    theta_up[i,]<-theta[i,,laststep[i]] 
    setTxtProgressBar(pb, i)
  }
  return(list(Indicator=Indicator,theta_up=theta_up,laststep=laststep)) 
} 

#####################################
psi<-function(x,tau)
{
  tau-1*(x<0)
}
rk<-function(T,k,ehat0,ehat,tau,weight)
{ 
  #tau is a scalar; k is the lag
  #ehat0 are residuals from estimation; ehat are residuals or bootstrap residuals   
  #weight: (length(ehat)-k)*1 vector of equal weight in diagnosis or the random weight in bootstrap 
  var_a=var(abs(ehat0))
  rk=sum(weight*psi(ehat[(k+1):T],tau)*abs(ehat[1:(T-k)]))/(T)/sqrt((tau-tau^2)*var_a)
  
  return(rk)
}
CICQE<-function(Y)
{
  M<-100
  NN<-900
  a<-1
  b<-1000
  N<-30
  tau1<-0.1
  tau2<-0.05
  p<-1
  q<-1
  B<-1000
  #######################################################################
  # This package is used to forecast conditional quantiles for the data #         
  # and to construct their confidence intervals (CIs)                   # 
  ##############################Input####################################
  # tau is the quantile level of interest                               #
  # p and q are the orders in GARCH(p,q) model                          #
  # B is size of Bootstrap                                              #
  # a and b determine the starting and ending points respectively       #
  # M is the out-sample size of forecasting                             #
  # data is the t.s. data to construct CIs for conditional quantiles    #
  # NN=b-M #start point of prediction in original data                  #
  #######################################################################
  #########################################################################    
  cat("\nsize of rolling forecasts (total number of iterations) = 1635:\n")
  #matrices of results
  
  thetaH1<-thetaH2<-array(0,dim=c(N,p+q+1,M)) #hybrid and QMLE
  QuantileH1<-QuantileH2<-matrix(0,N,M) #one-step ahead CQE and its CI 
  #one-step ahead CQE in bootstrap
  QHY1<-QHY2<-matrix(0,N,M)
  
  NN=b-M  #start point of prediction in original data
  n0=max(p,q)
  
  Ind<-matrix(0,N,M)
  
  for (l in 1:M)
  {
    Xlog=Y[,a:(NN+l-1)]
    T=dim(Xlog)[2]
    X<-Xlog
    Theta_int<-matrix(0,N,p+q+1)
    
    for (i in 1:N) {
      
      ##############################Data generation############################  
      #sample size
      ##############################Estimation#################################
      ########################GARCH(1,1) hybrid estimation#####################
      ####step I. QMLE for GARCH(1,1) to obtain h_t
      spec1_1=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(p, q)),                 
                         mean.model = list(armaOrder = c(0, 0),include.mean = FALSE),                 
                         distribution.model = "norm") 
      XG<-X
      garch1_1=ugarchfit(spec = spec1_1, data = XG[i,], solver = 'hybrid', fit.control = list(stationarity = 1))
      Theta_int[i,]=as.vector(coef(garch1_1)) #initial value for QMLE
      intc=rowMeans(X[,1:5]^2) #initial value of X_t^2 and h_t for t<=0 
      #initial value can also be mean(X^2) or X[1]^2 or Theta_int[1]/(1-Theta_int[2]-Theta_int[3]) 
      #the influence of intial value is small
    }
    QMLE=QMLE_NR(N,T,MM=20,p,q,tolerance=10^(-6),X,Theta_int,intc)    
    QMLEInd=QMLE$Indicator
    thetatilde<-matrix(0,N,3)
    
    XX=array(0,dim=c(N,T,2))
    YY<-matrix(0,N,T)
    
    for (i in 1:N) {
      
      if(QMLEInd[i]!=0)
      {          
        Ind[i,l]=1#;print(Ind)
        thetatilde[i,]=Theta_int[i,]  #(p+q+1)*1 vector
      }else
      {
        Ind[i,l]=0#;print(Ind)
        thetatilde[i,]=QMLE$theta_up[i,]  #(p+q+1)*1 vector
      }    
      httilde=h_update(N,T,X,thetatilde,intc)  #length(httilde)=N  
      
      ####step II. WCQE at specific tau
      #Create Responose and Regressor 
      X1=cbind(intc,X[,1:(T-1)]^2)
      X2=cbind(intc,httilde[,1:(T-1)])
      #Regressor without intercept in step II of GARCH(1,1) model; n*(p+q)  
      
      XX[i,,]=cbind(X1[i,],X2[i,])
      YY[i,]=X[i,]^2*sign(X[i,]) #Responose
      fit1=rq(YY[i,] ~ XX[i,,], tau=tau1, weights=1/httilde[i,])
      fit2=rq(YY[i,] ~ XX[i,,], tau=tau2, weights=1/httilde[i,])
      
      thetaH1[i,,l]=as.vector(coef(fit1))
      thetaH2[i,,l]=as.vector(coef(fit2))
      
      #############################One-step ahead CQE############################
      QHY1[i,l]=as.numeric(crossprod(thetaH1[i,,l],c(1,X[i,T]^2,httilde[i,T])))
      QuantileH1[i,l]=sqrt(abs(QHY1[i,l]))*sign(QHY1[i,l])
      QHY2[i,l]=as.numeric(crossprod(thetaH2[i,,l],c(1,X[i,T]^2,httilde[i,T])))
      QuantileH2[i,l]=sqrt(abs(QHY2[i,l]))*sign(QHY2[i,l])
    }
  }
  return(list(QuantileH1=QuantileH1,QuantileH2=QuantileH2))
  #############################Output#######################################
  # QuantileH is the rolling forecast of the conditional quantile          #
  # CIL and CIU are lower and upper bounds for rolling forecasts           #
  ##########################################################################
}
###########################################
################################
library(parallel)
ncores<-100
index_chunks<-split(1:200,cut(1:200,ncores,labels = FALSE))
cl<-makeCluster(ncores)
clusterEvalQ(cl,{library(quantreg)
  library(rugarch)})
clusterExport(cl,c("CICQE", "h_update", "derivativeG", "QMLE_NR","psi","rk", "x"))
results001<-parLapply(cl,index_chunks,function(idx){lapply(idx,function(i)CICQE(x[,,i]))})
stopCluster(cl)
###############################################

results00101<-array(0,dim = c(30,100,200))
for (i in 1:100) {
  for (j in 1:2) {
    results00101[,,2*(i-1)+j]<-results001[[i]][[j]]$QuantileH1
  }
}
##############################################
results001005<-array(0,dim = c(30,100,200))
for (i in 1:100) {
  for (j in 1:2) {
    results001005[,,2*(i-1)+j]<-results001[[i]][[j]]$QuantileH2
  }
}
################################################
load("raw_results/Table1/results00101.RData")
load("raw_results/Table1/results001005.RData")
mae1<-mae2<-rep(0,30)
for (i in 1:30) {
  mae1[i]<-mean(abs(results001005[i,,]-truequantile005[i,901:1000,]))
}
for (i in 1:30) {
  mae2[i]<-mean(abs(results00101[i,,]-truequantile01[i,901:1000,]))
}
mae1
mae2
#########################################################
