return<-read.csv("final return 2016.csv",header = TRUE)
return<-return[,-1]
return0<-ts(return)
dim(return0)
return0<-t(return0)
N<-71

library(BigVAR)
library(frequencyConnectedness)
big_var_est<- function(Y) {Model1 =constructModel(as.matrix(Y),1,struct="Basic",gran=c(50,10),
                                                  h = 20,cv = "Rolling",verbose = FALSE,IC = TRUE,VARX = list(), 
                                                  T1 = floor(nrow(Y)/3),T2 = floor(2 * nrow(Y)/3),ONESE = FALSE,
                                                  ownlambdas = FALSE,recursive = FALSE,dates = as.character(NULL), 
                                                  window.size = 500,separate_lambdas = FALSE)
MYResults= cv.BigVAR(Model1)
return(MYResults)}
allw<-array(NA,dim = c(N,N,294))
for (i in 1:294) {
  wdata<-return[1:(1454+i),]
  ww<-big_var_est(wdata)
  dynamic<-spilloverDY12(ww,n.ahead = 100,no.corr = F)
  allw[,,i]<-dynamic$tables[[1]]}

X<-return0
M<-294
NN<-1455
a<-1
b<-1749
p<-1
q<-1
B<-1000
Y<-return0
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
#####################################

CICQEmean<-function(N,allw,p=1,q=1,a,b,M,Y)
{
 
  ##############################Input####################################
  # p and q are the orders in GARCH(p,q) model                          #
  # B is size of Bootstrap                                              #
  # a and b determine the starting and ending points respectively       #
  # M is the out-sample size of forecasting                             #
  # data is the t.s. data to construct CIs for conditional quantiles    #
  # NN=b-M #start point of prediction in original data                  #
  #######################################################################
  #########################################################################    

  
  thetaH<-array(0,dim=c(N,p+q+1+1,M)) #hybrid and QMLE
  QuantileH<-matrix(0,N,M) #one-step ahead CQE and its CI 
  #one-step ahead CQE in bootstrap
  QHY<-matrix(0,N,M)
  
  NN=b-M  #start point of prediction in original data
  n0=max(p,q)
  Ind<-matrix(0,N,M)
  
  for (l in 1:M)
  {
    Xlog=Y[,a:(NN+l-1)]
    T=dim(Xlog)[2]
    X<-Xlog
    W<-allw[,,l]
    Theta_int<-matrix(0,N,p+q+1+1)
    
    for (i in 1:N) {
      
      ##############################Data generation############################  
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
    }
    httilde=h_update(N,T,W,X,thetatilde,intc) 
    ####step II. WCQE at specific tau
    #Create Responose and Regressor 
    X1=cbind(intc,X[,1:(T-1)]^2)
    X2=cbind(intc,httilde[,1:(T-1)])
    X3=W%*%X1
    #Regressor without intercept in step II of GARCH(1,1) model; n*(p+q)  
    for (i in 1:N) {
      thetaH[i,,l]=thetatilde[i,]
      #############################One-step ahead CQE############################
      QHY[i,l]=as.numeric(crossprod(thetaH[i,,l],c(1,X[i,T]^2,httilde[i,T],W[i,]%*%(X[,T]^2))))
      QuantileH[i,l]=sqrt(abs(QHY[i,l]))*sign(QHY[i,l])
    }
  }
  return(list(QuantileH=QuantileH))
  #############################Output#######################################

  ##########################################################################
}

run_parallel_cicqe <- function(task) {
  # 子进程环境初始化
  library(quantreg)
  library(rugarch)
  
  # 执行计算
  # 注意：参数直接从 task 列表中提取
  res <- CICQEmean(
    N = 71, 
    allw = task$allw,
    p = 1, q = 1, 
    a = 1, 
    b = task$b, 
    M = 21, 
    Y = task$Y
  )
  return(res)
}
library(parallel)
cl<-makeCluster(14)
clusterEvalQ(cl,{
  library(quantreg)
  library(rugarch)})
tasks <- list(
  task1 = list(allw=allw[,,1:21], b = 1476, Y = return0[, 1:1476]),
  task2 = list(allw=allw[,,22:42], b = 1497, Y = return0[, 1:1497]),
  task3 = list(allw=allw[,,43:63], b = 1518, Y = return0[, 1:1518]),
  task4 = list(allw=allw[,,64:84], b = 1539, Y = return0[, 1:1539]),
  task5 = list(allw=allw[,,85:105], b = 1560, Y = return0[, 1:1560]),
  task6 = list(allw=allw[,,106:126], b = 1581, Y = return0[, 1:1581]),
  task7 = list(allw=allw[,,127:147], b = 1602, Y = return0[, 1:1602]),
  task8 = list(allw=allw[,,148:168], b = 1623, Y = return0[, 1:1623]),
  task9 = list(allw=allw[,,169:189], b = 1644, Y = return0[, 1:1644]),
  task10 = list(allw=allw[,,190:210], b = 1665, Y = return0[, 1:1665]),
  task11 = list(allw=allw[,,211:231], b = 1686, Y = return0[, 1:1686]),
  task12 = list(allw=allw[,,232:252], b = 1707, Y = return0[, 1:1707]),
  task13 = list(allw=allw[,,253:273], b = 1728, Y = return0[, 1:1728]),
  task14 = list(allw=allw[,,274:294], b = 1749, Y = return0[, 1:1749])
)

clusterExport(cl,varlist = c("CICQEmean", "h_update", "h_update1", "derivativeG", "QMLE_NR", 
                             "psi", "rk", "allw", "return0"))
results<-clusterApplyLB(cl,tasks,run_parallel_cicqe)              
stopCluster(cl) 
nqforecastmean<-cbind(results[[1]]$QuantileH,results[[2]]$QuantileH,results[[3]]$QuantileH,
                    results[[4]]$QuantileH,results[[5]]$QuantileH,results[[6]]$QuantileH,
                    results[[7]]$QuantileH,results[[8]]$QuantileH,results[[9]]$QuantileH,
                    results[[10]]$QuantileH,results[[11]]$QuantileH,results[[12]]$QuantileH,
                    results[[13]]$QuantileH,results[[14]]$QuantileH)
write.csv(nqforecastmean,"nqgarch-forecastmean.csv")

###############################




###############################################
CICQE<-function(N,allw,tau,p=1,q=1,B,a,b,M,Y)
{
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
  
  thetaH<-array(0,dim=c(N,p+q+1+1,M)) #hybrid and QMLE
  QuantileH<-matrix(0,N,M) #one-step ahead CQE and its CI 
  #one-step ahead CQE in bootstrap
  QHY<-matrix(0,N,M)
  
  NN=b-M  #start point of prediction in original data
  n0=max(p,q)
  Ind<-matrix(0,N,M)
  
  for (l in 1:M)
  {
    Xlog=Y[,a:(NN+l-1)]
    T=dim(Xlog)[2]
    X<-Xlog
    W<-allw[,,l]
    Theta_int<-matrix(0,N,p+q+1+1)
    
    for (i in 1:N) {
      
      ##############################Data generation############################  
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
    }
    httilde=h_update(N,T,W,X,thetatilde,intc)  #length(httilde)=N  
    
    ####step II. WCQE at specific tau
    
    #Create Responose and Regressor 
    X1=cbind(intc,X[,1:(T-1)]^2)
    X2=cbind(intc,httilde[,1:(T-1)])
    X3=W%*%X1
    #Regressor without intercept in step II of GARCH(1,1) model; n*(p+q)  
    for (i in 1:N) {
      XX[i,,]=cbind(X1[i,],X2[i,],X3[i,])
      YY[i,]=X[i,]^2*sign(X[i,]) #Responose
      fit=rq(YY[i,] ~ XX[i,,], tau=tau, weights=1/httilde[i,])
      thetaH[i,,l]=as.vector(coef(fit))   
      #############################One-step ahead CQE############################
      QHY[i,l]=as.numeric(crossprod(thetaH[i,,l],c(1,X[i,T]^2,httilde[i,T],W[i,]%*%(X[,T]^2))))
      QuantileH[i,l]=sqrt(abs(QHY[i,l]))*sign(QHY[i,l])
    }
    
  }
  return(list(QuantileH=QuantileH))
  #############################Output#######################################
  # QuantileH is the rolling forecast of the conditional quantile          #
  ##########################################################################
}
run_parallel_cicqe <- function(task) {
  # 子进程环境初始化
  library(quantreg)
  library(rugarch)
  
  # 执行计算
  # 注意：参数直接从 task 列表中提取
  res <- CICQE(
    N = 71, 
    allw = task$allw,
    tau = 0.5, 
    p = 1, q = 1, 
    B = 1000, 
    a = 1, 
    b = task$b, 
    M = 21, 
    Y = task$Y
  )
  return(res)
}
library(parallel)
cl<-makeCluster(14)
clusterEvalQ(cl,{
  library(quantreg)
  library(rugarch)})
tasks <- list(
  task1 = list(allw=allw[,,1:21], b = 1476, Y = return0[, 1:1476]),
  task2 = list(allw=allw[,,22:42], b = 1497, Y = return0[, 1:1497]),
  task3 = list(allw=allw[,,43:63], b = 1518, Y = return0[, 1:1518]),
  task4 = list(allw=allw[,,64:84], b = 1539, Y = return0[, 1:1539]),
  task5 = list(allw=allw[,,85:105], b = 1560, Y = return0[, 1:1560]),
  task6 = list(allw=allw[,,106:126], b = 1581, Y = return0[, 1:1581]),
  task7 = list(allw=allw[,,127:147], b = 1602, Y = return0[, 1:1602]),
  task8 = list(allw=allw[,,148:168], b = 1623, Y = return0[, 1:1623]),
  task9 = list(allw=allw[,,169:189], b = 1644, Y = return0[, 1:1644]),
  task10 = list(allw=allw[,,190:210], b = 1665, Y = return0[, 1:1665]),
  task11 = list(allw=allw[,,211:231], b = 1686, Y = return0[, 1:1686]),
  task12 = list(allw=allw[,,232:252], b = 1707, Y = return0[, 1:1707]),
  task13 = list(allw=allw[,,253:273], b = 1728, Y = return0[, 1:1728]),
  task14 = list(allw=allw[,,274:294], b = 1749, Y = return0[, 1:1749])
)

clusterExport(cl,varlist = c("CICQE", "h_update", "h_update1", "derivativeG", "QMLE_NR", 
                             "psi", "rk", "allw", "return0"))
results<-clusterApplyLB(cl,tasks,run_parallel_cicqe)              
stopCluster(cl) 
nq05forecast<-cbind(results[[1]]$QuantileH,results[[2]]$QuantileH,results[[3]]$QuantileH,
                      results[[4]]$QuantileH,results[[5]]$QuantileH,results[[6]]$QuantileH,
                      results[[7]]$QuantileH,results[[8]]$QuantileH,results[[9]]$QuantileH,
                      results[[10]]$QuantileH,results[[11]]$QuantileH,results[[12]]$QuantileH,
                      results[[13]]$QuantileH,results[[14]]$QuantileH)
write.csv(nq05forecast,"nqgarch-forecast05.csv")
###################################################################
###############



###############################
forecast05<-read.csv(file = "nqgarch-forecast05.csv",header = TRUE)
forecastmean<-read.csv(file = "nqgarch-forecastmean.csv", header = TRUE)
observation<-return0[,1456:1749]
forecast05<-forecast05[,-1]
forecastmean<-forecastmean[,-1]
e05<-forecast05-observation
emean<-forecastmean-observation
library(forecast)
setest<-matrix(0,71,2)
for (i in 1:71) {
  dm<-dm.test(as.numeric(e05[i,]),as.numeric(emean[i,]),alternative="less",power=2)
  setest[i,1]<-dm$statistic
  setest[i,2]<-dm$p.value  
}
round(setest,4)
aetest<-matrix(0,71,2)
for (i in 1:71) {
  dm<-dm.test(as.numeric(e05[i,]),as.numeric(emean[i,]),alternative="less",power=1)
  aetest[i,1]<-dm$statistic
  aetest[i,2]<-dm$p.value  
}
round(aetest,4)

################################################################
e05<-ts(e05)
library(moments)
library(tseries)
e05stat<-matrix(0,71,3)
for (i in 1:71) {
  e05stat[i,1]<-skewness(e05[i,])
  e05stat[i,2]<-kurtosis(e05[i,])
  e05stat[i,3]<-jarque.bera.test(e05[i,])$p.value
}
round(e05stat,4)
emean<-ts(emean)
emeanstat<-matrix(0,71,3)
for (i in 1:71) {
  emeanstat[i,1]<-skewness(emean[i,])
  emeanstat[i,2]<-kurtosis(emean[i,])
  emeanstat[i,3]<-jarque.bera.test(emean[i,])$p.value
}
round(emeanstat,4)
#############################################################################

