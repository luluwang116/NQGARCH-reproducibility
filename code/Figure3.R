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
######################

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
CQEst<-function(tau,p=1,q=1,X,W)
{
  #######################################################################
  # This function is used to perform estimation based on the hybrid     #
  # estimation procedure and do diagnostic checking for conditional     #
  # quantiles based on the QACF of residuals                            #
  #############################Input#####################################
  # tau is the quantile level of interest                               #
  # p and q are the orders in GARCH(p,q) model                          #
  # L is a predetermined integer for the maximum lag in QACF            #
  # X is the t.s. data for analysis                                     #
  # B is size of Bootstrap                                              #
  #######################################################################
  N=nrow(X) #N is the spatial sample size
  T=ncol(X) #N is the time series length
  ##############################Matrices of results#######################
  
  ##############################QMLE for GARCH(1,1)##########################
  Theta_int<-matrix(NA,N,4)
  for (i in 1:N) {
    spatial_lag<-t(W%*%(X^2))[-T,]
    spec1_1=ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1),external.regressors=as.matrix(spatial_lag[,i])),                 
                       mean.model = list(armaOrder = c(0, 0),include.mean = FALSE),                 
                       distribution.model = "norm")
    ts<-ts(X[i,-1])
    garch1_1=ugarchfit(spec = spec1_1, data =ts , solver = 'hybrid', fit.control = list(stationarity = 1))
    Theta_int[i,]<-as.vector(coef(garch1_1)) 
  }
  intc=rowMeans(X[,1:5]^2) #initial value of X_t^2 and h_t for t<=0 
  #initial value can also be mean(X^2) or X[1]^2 or Theta_int[1]/(1-Theta_int[2]-Theta_int[3]) 
  #the influence of intial value is small
  bad <- which(
    Theta_int[,2] < 1e-5 |
      Theta_int[,3] > 0.95 |
      Theta_int[,4] < 1e-5 |
      rowSums(Theta_int[,2:4]) >= 0.95
  )
  
  if (length(bad) > 0) {
    target <- apply(
      Theta_int[-bad,2:4,drop=FALSE],
      2, median )
    Theta_int[bad,2:4] <-
      0.5*Theta_int[bad,2:4,drop=FALSE] +
      0.5*matrix(target,nrow=length(bad),ncol=3,byrow=TRUE)
    Theta_int[bad,1] <- pmax(( 1-rowSums(Theta_int[bad,2:4,drop=FALSE]))*intc[bad],1e-5) }
  
  QMLE=QMLE_NR(N,T,W,MM=20,p,q,tolerance=10^(-6),X,Theta_int,intc)    
  QMLEInd=QMLE$Indicator 
  thetatilde<-matrix(0,N,4)
  Ind<-rep(0,N)
  Y<-matrix(0,N,T)
  for (i in 1:N) {
    
    if(QMLEInd[i]!=0)
    {          
      Ind[i]=1#;print(Ind)
      thetatilde[i,]=Theta_int[i,]  #(p+q+1+1)*1 vector
    }else
    {
      Ind[i]=0#;print(Ind)
      thetatilde[i,]=QMLE$theta_up[i,]  #(p+q+1)*1 vector
    }
    ######################GARCH(p,q) two-step estimation#######################  
    Y[i,]=X[i,]^2*sign(X[i,]) #Responose
  }
  #####step I. QMLE for GARCH(1,1) to obtain h_t
  httilde=h_update(N,T,W,X,thetatilde,intc)  #length(httilde)=N  
  #####step II. WCQE at specific tau
  
  #Create Regressor
  X1=cbind(intc,X[,1:(T-1)]^2)
  X2=cbind(intc,httilde[,1:(T-1)])
  X3=W%*%X1
  XX=array(0,dim=c(N,T,3)) #Regressor without intercept in step II of GARCH(1,1) model; N*(p+q)  
  
  
  thetahat<-matrix(0,N,4)
  Der=derivativeG(N,T,W,X,thetatilde,intc)
  P1h=Der$P1h #N*T*(p+q+1)
  P1ell=Der$P1ell #N*T*(p+q+1)
  colnames(thetahat)=c("alpha0tau","alpha1tau","beta1tau","rhotau")
  
  pb <- txtProgressBar(min = 0, max = N,style = 3, width = 50, char = "=")
  for (i in 1:N) {
    
    
    XX[i,,]=cbind(X1[i,],X2[i,],X3[i,])
    fit=rq(Y[i,] ~ XX[i,,], tau=tau,weights=1/httilde[i,])
    thetahat[i,]=as.vector(coef(fit))
  }
  return(list(thetatilde=thetatilde,thetahat=thetahat)
  )
  
  #############################Output#######################################
  # thetatilde and sdQMLE are estimates and standard errors of QMLE        #
  # thetahat and sdE are estimates and standard errors of hybrid estimator #
  
  ##########################################################################   
}

##########################################################################   
##########################################################################   
set.seed(212745678)##tau=0.05, T=1000 and student t5
n <- 30
T <- 1050
S <- 500

xit <- array(0, dim = c(n, T, S))
eta <- array(0, dim = c(n, T, S))
h   <- array(0, dim = c(n, T, S))

############################################################

eta[] <- rt(
  length(eta),df=5)/(sqrt(5/3))

############################################################
# Generate Network-GARCH process
############################################################

for(s in 1:S){
  
  h[, 1, s] <-
    alpha0 +
    alpha1 * h0 +
    beta * h0 +
    rho * as.numeric(
      W %*% h0
    )
  
  xit[, 1, s] <-
    sqrt(h[, 1, s]) *
    eta[, 1, s]
  
  
  for(t in 2:T){
    
    x2_lag <- xit[, t - 1, s]^2
    
    h[, t, s] <-
      alpha0 +
      alpha1 * x2_lag +
      beta * h[, t - 1, s] +
      rho * as.numeric(
        W %*% x2_lag
      )
    
    xit[, t, s] <-
      sqrt(h[, t, s]) *
      eta[, t, s]
  }
}


############################################################
# Remove first 50 observations
############################################################

xit <- xit[, 51:T, , drop = FALSE]


dim(xit)

######################
library(parallel)

ncores <- 100
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(quantreg)
  library(rugarch)
})

clusterExport(
  cl,
  varlist = c("CQEst","h_update", "h_update1", "derivativeG", "QMLE_NR", 
              "psi", "rk", "xit", "W"),
  envir = environment()
)

results_coef <- parLapplyLB(
  cl,
  1:500,
  function(s) {
    fit <- CQEst(
      tau = 0.05,
      p = 1,
      q = 1,
      X = xit[,,s],
      W = W
    )
    
    fit$thetahat
  }
)

stopCluster(cl)
############################################
coef <- simplify2array(results_coef)

dim(coef)
garchcoeftau<--(qt(0.05,df=5)/sqrt(5/3))^2*garchcoef
bias1<-matrix(0,30,4)
for (i in 1:30) {
  bias1[i,]<-rowMeans(coef[i,,])-garchcoeftau[i,]
}
bias1
############################################################
set.seed(212745679)##tau=0.05, T=1500 and student t5
n <- 30
T <- 1550
S <- 500

xit <- array(0, dim = c(n, T, S))
eta <- array(0, dim = c(n, T, S))
h   <- array(0, dim = c(n, T, S))

############################################################

eta[] <- rt(
  length(eta),df=5)/(sqrt(5/3))

############################################################
# Generate Network-GARCH process
############################################################

for(s in 1:S){
  
  h[, 1, s] <-
    alpha0 +
    alpha1 * h0 +
    beta * h0 +
    rho * as.numeric(
      W %*% h0
    )
  
  xit[, 1, s] <-
    sqrt(h[, 1, s]) *
    eta[, 1, s]
  
  
  for(t in 2:T){
    
    x2_lag <- xit[, t - 1, s]^2
    
    h[, t, s] <-
      alpha0 +
      alpha1 * x2_lag +
      beta * h[, t - 1, s] +
      rho * as.numeric(
        W %*% x2_lag
      )
    
    xit[, t, s] <-
      sqrt(h[, t, s]) *
      eta[, t, s]
  }
}


############################################################
# Remove first 50 observations
############################################################

xit <- xit[, 51:T, , drop = FALSE]


dim(xit)
######################
library(parallel)

ncores <- 100
cl <- makeCluster(ncores)

clusterEvalQ(cl, {
  library(quantreg)
  library(rugarch)
})

clusterExport(
  cl,
  varlist = c("CQEst","h_update", "h_update1", "derivativeG", "QMLE_NR", 
              "psi", "rk", "xit", "W"),
  envir = environment()
)

results_coef <- parLapplyLB(
  cl,
  1:500,
  function(s) {
    fit <- CQEst(
      tau = 0.05,
      p = 1,
      q = 1,
      X = xit[,,s],
      W = W
    )
    
    fit$thetahat
  }
)

stopCluster(cl)

###################################################
coef <- simplify2array(results_coef)

dim(coef)

garchcoeftau<--(qt(0.05,df=5)/sqrt(5/3))^2*garchcoef
bias2<-matrix(0,30,4)
for (i in 1:30) {
  bias2[i,]<-rowMeans(coef[i,,])-garchcoeftau[i,]
}
bias2
####################################################
bias0051000<-bias1
bias0051500<-bias2
alpha0005<-cbind(bias0051000[,1],bias0051500[,1])
colnames(alpha0005)<-c("T=1000","T=1500")
boxplot(alpha0005,col="skyblue",xlab="alpha0",ylab="Bias",medcol="royalblue",border="skyblue")
alpha1005<-cbind(bias0051000[,2],bias0051500[,2])
colnames(alpha1005)<-c("T=1000","T=1500")
boxplot(alpha1005,col="skyblue",xlab="alpha1",ylab="Bias",medcol="royalblue",border="skyblue")
beta005<-cbind(bias0051000[,3],bias0051500[,3])
colnames(beta005)<-c("T=1000","T=1500")
boxplot(beta005,col="skyblue",xlab="beta",ylab="Bias",medcol="royalblue",border="skyblue")
lambda005<-cbind(bias0051000[,4],bias0051500[,4])
colnames(lambda005)<-c("T=1000","T=1500")
boxplot(lambda005,col="skyblue",xlab="lambda",ylab="Bias",medcol="royalblue",border="skyblue")
############################################