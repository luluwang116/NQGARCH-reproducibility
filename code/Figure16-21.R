Y <-read.csv("final return 2016.csv")
Y<-as.matrix(Y[1:818,2:72])
library(BigVAR)
library(frequencyConnectedness)
big_var_est<- function(Y) {Model1 = constructModel(as.matrix(Y),1,struct="Basic",gran=c(50,10),
                                                   h = 20,cv = "Rolling",verbose = FALSE,IC = TRUE,VARX = list(), 
                                                   T1 = floor(nrow(Y)/3),T2 = floor(2 * nrow(Y)/3),ONESE = FALSE,
                                                   ownlambdas = FALSE,recursive = FALSE,dates = as.character(NULL), 
                                                   window.size = 500,separate_lambdas = FALSE)
MYResults= cv.BigVAR(Model1)}
oo <- big_var_est(Y)
spilloverDY12(oo, n.ahead = 100, no.corr = F)
staticnet=spilloverDY12(oo, n.ahead = 100, no.corr = F)
W1<-as.matrix(staticnet$tables[[1]])
write.csv(W1,"subsamplematrix1.csv")
############################################
Y <-read.csv("final return 2016.csv")
Y<-as.matrix(Y[1177:1749,2:72])
library(BigVAR)
library(frequencyConnectedness)
big_var_est<- function(Y) {Model1 = constructModel(as.matrix(Y),1,struct="Basic",gran=c(50,10),
                                                   h = 20,cv = "Rolling",verbose = FALSE,IC = TRUE,VARX = list(), 
                                                   T1 = floor(nrow(Y)/3),T2 = floor(2 * nrow(Y)/3),ONESE = FALSE,
                                                   ownlambdas = FALSE,recursive = FALSE,dates = as.character(NULL), 
                                                   window.size = 500,separate_lambdas = FALSE)
MYResults= cv.BigVAR(Model1)}
oo <- big_var_est(Y)
spilloverDY12(oo, n.ahead = 100, no.corr = F)
staticnet=spilloverDY12(oo, n.ahead = 100, no.corr = F)
W2<-as.matrix(staticnet$tables[[1]])
write.csv(W2,"subsamplematrix2.csv")
cor(as.vector(W1), as.vector(W2))
########################################
return<-read.csv("final return 2016.csv",header = TRUE)
return<-return[,-1]
return0<-ts(return)
dim(return0)
return0<-t(return0)

N<-71
T1<-818
T2<-573
p<-1
q<-1
B<-1000
#################################
library(rugarch)
library(quantreg)
#########################

##############################
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
CQEst<-function(tau,p=1,q=1,L,X,B,W)
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
  Sigma1hatB<-array(0,dim=c(N,p+q+1+1,p+q+1+1))  #bootstrapped covariance matrices for WCQE
  cvQ<-matrix(0,N,L)#critical values for Q
  
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
  ehat<-matrix(0,N,T)
  R<-matrix(0,N,L)
  Jtilde<-array(0,dim=c(N,4,4))
  Jinv<-array(0,dim=c(N,4,4))
  
  WW<-array(0,dim=c(N,T,B))
  e_c<-matrix(0,N,4)
  thetatildeB<-matrix(0,N,p+q+1+1)  #QMLE after correction in bootstrapping
  httildeB<-matrix(0,N,T)
  thetahatB<-matrix(0,N,p+q+1+1)  #WCQE in bootstrapping
  dE<-array(0,dim=c(N,p+q+1+1,B)) #sqrt(T)*(thetahatB-thetahat)
  ehatB<-matrix(0,N,T)
  RB<-matrix(0,N,L) #L*1 vector for bootstrap test statistics
  dT<-array(0,dim=c(N,L,B)) #sqrt(N)*(RhatB-Rhat)
  
  X2B<-matrix(0,N,T)
  XXB=array(0,dim=c(N,T,3))
  
  kappa<-matrix(0,N,1)
  sdthetatilde<-matrix(0,N,4)
  sdthetahatB<-matrix(0,N,p+q+1+1)  #bootstrapped standard errors for WCQE
  Sigma2hatB<-array(0,dim=c(N,L,L))  #bootstrapped covariance matrices for R
  sdRB<-matrix(0,N,L)  #bootstrapped standard errors for R
  LCIRB<-UCIRB<-matrix(0,nrow=N,ncol=L)  #confidence intervals for R
  QR<-matrix(0,N,1)
  
  
  
  
  
  set.seed(1234567)
  
  Der=derivativeG(N,T,W,X,thetatilde,intc)
  P1h=Der$P1h #N*T*(p+q+1)
  P1ell=Der$P1ell #N*T*(p+q+1)
  colnames(thetahat)=c("alpha0tau","alpha1tau","beta1tau","rhotau")
  
  pb <- txtProgressBar(min = 0, max = N,style = 3, width = 50, char = "=")
  for (i in 1:N) {
    
    
    XX[i,,]=cbind(X1[i,],X2[i,],X3[i,])
    fit=rq(Y[i,] ~ XX[i,,], tau=tau, weights=1/httilde[i,])
    thetahat[i,]=as.vector(coef(fit))  
    ######################GARCH(p,q) diagnostic checking#######################
    ehat[i,]=(Y[i,]-t(as.matrix(cbind(1,XX[i,,])%*%thetahat[i,])))/httilde[i,]
    for (k in 1:L)
    {
      R[i,k]=rk(T,k,ehat[i,],ehat[i,],tau,weight=rep(1,T-k)) 
    }
    #######################Correction matrices##################################
    
    Jtilde[i,,]=crossprod(P1h[i,,]/httilde[i,])/T  #(p+q+1)*(p+q+1)
    Jinv[i,,]=solve(Jtilde[i,,], tol = 1e-16)
    ##################################Bootstrap#################################
    
  }
  for (i in 1:N) {
    
    
    WW[i,,]<-matrix(rexp(T*B, rate = 1),T,B)  # standard exponential
    
    b=1
    while (b<=B)
    {
      ##########################Weights generation###############################
      weight=WW[i,,b]  #standard exponential
      ############################Correction terms###############################
      e_c[i,]=Jinv[i,,]%*%apply((weight-1)*P1ell[i,,],2,mean)  #correction term
      ####step I. corrected QMLE for GARCH(1,1) to obtain h_tB for bootstrap
      thetatildeB[i,]=thetatilde[i,]-e_c[i,] #corrected QMLE  
      
      httildeB[i,]=h_update1(T,W[i,],X[i,],X,thetatildeB[i,],intc[i],intc)   
      ####step II. WCQE at specific tau for bootstrap    
      X2B[i,]=c(intc[i],httildeB[i,1:(T-1)])
      XXB[i,,]=cbind(X1[i,],X2B[i,],X3[i,]) #Regressor without intercept in step II for bootstrap  
      
      fitB=rq(Y[i,] ~ XXB[i,,], tau=tau, weights=weight/httilde[i,])
      thetahatB[i,]=as.vector(coef(fitB))       
      dE[i,,b]=sqrt(T)*(thetahatB[i,]-thetahat[i,])   
      ##############GARCH(p,q) diagnostic checking for bootstrap##############
      ehatB[i,]=(Y[i,]-t(as.matrix(cbind(1,XXB[i,,])%*%thetahatB[i,])))/httilde[i,]
      
      for (k in 1:L)
      {
        RB[i,k]=rk(T,k,ehat[i,],ehatB[i,],tau,weight=weight[(k+1):T])
      }    
      dT[i,,b]=sqrt(T)*(RB[i,]-R[i,])
      
      b=b+1
      #print(b) 
    } 
    
    #QMLE
    kappa[i,]=mean(X[i,]^4/httilde[i,]^2)
    sdthetatilde[i,]=sqrt(diag((kappa[i,]-1)*Jinv[i,,]/T))
    
    #WCQE
    sdthetahatB[i,]=apply(dE[i,,],1,sd)/sqrt(T)
    
    #R
    Sigma2hatB[i,,]=var(t(dT[i,,]))
    sdRB[i,]=sqrt(diag(Sigma2hatB[i,,])/T) #apply(dT,1,sd)/sqrt(N)
    
    #check test statistic individually
    LCIRB[i,]=apply(dT[i,,],1,quantile,probs=0.025)/sqrt(T)  
    UCIRB[i,]=apply(dT[i,,],1,quantile,probs=0.975)/sqrt(T)
    
    #check test statistic jointly  
    QR[i,]=T*t(R[i,])%*%solve(Sigma2hatB[i,,])%*%(R[i,])
    setTxtProgressBar(pb, i)
  }
  return(list(thetatilde=thetatilde,sdQMLE=sdthetatilde,
              thetahat=thetahat,sdE=sdthetahatB,
              R=R,sdR=sdRB,QR=QR,            
              LCIRB=LCIRB,UCIRB=UCIRB)
  )
  
  #############################Output#######################################
  # thetatilde and sdQMLE are estimates and standard errors of QMLE        #
  # thetahat and sdE are estimates and standard errors of hybrid estimator #

  ##########################################################################   
}
###############################
####################################################
qfit1<-CQEst(tau=0.1,p=1,q=1,L=30,X=return0[,1:818],B=1000,W=W1)
phybrid1<-2*pt(-abs(qfit1$thetahat/qfit1$sdE),df=T1-4)
write.csv(qfit1$thetahat,file="coefrobust1-tau01.csv")
write.csv(phybrid1,file="pvaluerobust1-tau01.csv")
######################################################
qfit2<-CQEst(tau=0.25,p=1,q=1,L=30,X=return0[,1:818],B=1000,W=W1)
phybrid2<-2*pt(-abs(qfit2$thetahat/qfit2$sdE),df=T1-4)
write.csv(qfit2$thetahat,file="coefrobust1-tau025.csv")
write.csv(phybrid2,file="pvaluerobust1-tau025.csv")
######################################################
qfit3<-CQEst(tau=0.5,p=1,q=1,L=30,X=return0[,1:818],B=1000,W=W1)
phybrid3<-2*pt(-abs(qfit3$thetahat/qfit3$sdE),df=T1-4)
write.csv(qfit3$thetahat,file="coefrobust1-tau05.csv")
write.csv(phybrid3,file="pvaluerobust1-tau05.csv")
#######################################################
qfit4<-CQEst(tau=0.75,p=1,q=1,L=30,X=return0[,1:818],B=1000,W=W1)
phybrid4<-2*pt(-abs(qfit4$thetahat/qfit4$sdE),df=T1-4)
write.csv(qfit4$thetahat,file="coefrobust1-tau075.csv")
write.csv(phybrid4,file="pvaluerobust1-tau075.csv")
#######################################################
qfit5<-CQEst(tau=0.9 ,p=1,q=1,L=30,X=return0[,1:818],B=1000,W=W1)
phybrid5<-2*pt(-abs(qfit5$thetahat/qfit5$sdE),df=T1-4)
write.csv(qfit5$thetahat,file="coefrobust1-tau09.csv")
write.csv(phybrid5,file="pvaluerobust1-tau09.csv")
#######################################################
####################################################
qfit1<-CQEst(tau=0.1,p=1,q=1,L=30,X=return0[,1177:1749],B=1000,W=W2)
phybrid1<-2*pt(-abs(qfit1$thetahat/qfit1$sdE),df=T2-4)
write.csv(qfit1$thetahat,file="coefrobust2-tau01.csv")
write.csv(phybrid1,file="pvaluerobust2-tau01.csv")
######################################################
qfit2<-CQEst(tau=0.25,p=1,q=1,L=30,X=return0[,1177:1749],B=1000,W=W2)
phybrid2<-2*pt(-abs(qfit2$thetahat/qfit2$sdE),df=T2-4)
write.csv(qfit2$thetahat,file="coefrobust2-tau025.csv")
write.csv(phybrid2,file="pvaluerobust2-tau025.csv")
######################################################
qfit3<-CQEst(tau=0.5,p=1,q=1,L=30,X=return0[,1177:1749],B=1000,W=W2)
phybrid3<-2*pt(-abs(qfit3$thetahat/qfit3$sdE),df=T2-4)
write.csv(qfit3$thetahat,file="coefrobust2-tau05.csv")
write.csv(phybrid3,file="pvaluerobust2-tau05.csv")
#######################################################
qfit4<-CQEst(tau=0.75,p=1,q=1,L=30,X=return0[,1177:1749],B=1000,W=W2)
phybrid4<-2*pt(-abs(qfit4$thetahat/qfit4$sdE),df=T2-4)
write.csv(qfit4$thetahat,file="coefrobust2-tau075.csv")
write.csv(phybrid4,file="pvaluerobust2-tau075.csv")
#######################################################
qfit5<-CQEst(tau=0.9 ,p=1,q=1,L=30,X=return0[,1177:1749],B=1000,W=W2)
phybrid5<-2*pt(-abs(qfit5$thetahat/qfit5$sdE),df=T2-4)
write.csv(qfit5$thetahat,file="coefrobust2-tau09.csv")
write.csv(phybrid5,file="pvaluerobust2-tau09.csv")
#######################################################
coef01<-read.csv("coefrobust1-tau01.csv",header = TRUE)
coef025<-read.csv("coefrobust1-tau025.csv",header = TRUE)
coef05<-read.csv("coefrobust1-tau05.csv",header = TRUE)
coef075<-read.csv("coefrobust1-tau075.csv",header = TRUE)
coef09<-read.csv("coefrobust1-tau09.csv",header = TRUE)
pvalue01<-read.csv("pvaluerobust1-tau01.csv", header = TRUE)
pvalue025<-read.csv("pvaluerobust1-tau025.csv", header = TRUE)
pvalue05<-read.csv("pvaluerobust1-tau05.csv", header = TRUE)
pvalue075<-read.csv("pvaluerobust1-tau075.csv", header = TRUE)
pvalue09<-read.csv("pvaluerobust1-tau09.csv", header = TRUE)
coefalpha1<-cbind(coef01[,3],coef025[,3],coef05[,3],coef075[,3],coef09[,3])
coefbeta<-cbind(coef01[,4],coef025[,4],coef05[,4],coef075[,4],coef09[,4])
coeflambda<-cbind(coef01[,5],coef025[,5],coef05[,5],coef075[,5],coef09[,5])
pvaluealpha1<-cbind(pvalue01[,3],pvalue025[,3],pvalue05[,3],pvalue075[,3],pvalue09[,3])
pvaluebeta<-cbind(pvalue01[,4],pvalue025[,4],pvalue05[,4],pvalue075[,4],pvalue09[,4])
pvaluelambda<-cbind(pvalue01[,5],pvalue025[,5],pvalue05[,5],pvalue075[,5],pvalue09[,5])

alpha1<-matrix(NA,71,5)
for (i in 1:71) {
  for (j in 1:5) {
    if(pvaluealpha1[i,j]<0.1){
      alpha1[i,j]<-coefalpha1[i,j]
    }else{alpha1[i,j]<-0}
  }
}
beta<-matrix(NA,71,5)
for (i in 1:71) {
  for (j in 1:5) {
    if(pvaluebeta[i,j]<0.1){
      beta[i,j]<-coefbeta[i,j]
    }else{beta[i,j]<-0}
  }
}
lambda<-matrix(NA,71,5)
for (i in 1:71) {
  for (j in 1:5) {
    if(pvaluelambda[i,j]<0.1){
      lambda[i,j]<-coeflambda[i,j]
    }else{lambda[i,j]<-0}
  }
}


library(ggplot2)
library(dplyr)

commodities<-c("Corn, China", "Wheat, China", "Japonica Rice, China", "Pork, China", "Beef, China", 
               "Mutton, China", "Cotton, China", "Soybean, China", "Soybean Meal, China", 
               "Soybean Oil, China", "Rapeseed Meal, China", "Rapeseed Oil, China", "Palm Oil, China",
               "Pure Benzene, China", "Toluene, China", "O-Xylene, China", "Styrene, China", 
               "Butadiene, China", "P-Xylene, China", "PTA, China", "Epoxy Ethane, China", "Propylene, China", 
               "Polystyrene, China", "Natural Rubber, China", "PP, China", "Causticsoda, China", 
               "Potassic Fertilizers, China", "Xylene, China", "Glycol, China", "Acrylamide, China", 
               "Acrylic Acid, China", "Acrylonitrile, China", "Butyl Acrylate, China", "PET Chips, China", 
               "Methyl Methacrylate, China", "PVC, China", "Urea, China", "Methanol, China", 
               "Metallurgical Coke, China", "Hard Coking Coal, China", "Coking Coal, China", 
               "Hot-rolled Sheet/Coil, China", "Platinum, China", "Silver, China", "Palladium, China",
               "Gold, China", "Copper, China", "Nickel, China", "Aluminum, China", "Zinc, China",
               "Lead, China", "Tin, China", "Chrome Metal, China", "Silicon Metal, China", "Cobalt, China",
               "Ingot, China", "Electrolytic Manganese, China", "Nickel Sulfate, China", "Antimony, China",
               "Soybean Oil, US", "Diesel Oil, Singapore", "Gasoline, Japan", "Brent", "WTI", "Henry Hub Natural Gas",
               "Crude Oil, UAE Dubai", "Naphtha, Korea", "Gold, UK", "Silver, UK", "Platinum, UK", "Palladium, UK")
taus<-c("tau=0.1", "tau=0.25", "tau=0.5", "tau=0.75", "tau=0.9")
as.factor(taus)
data<-expand.grid(Commodity=commodities,Tau=taus)

dim(alpha1)<-355
data$alpha1<-alpha1
ggplot(data, aes(x=Tau, y=Commodity, fill=alpha1))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = c("red3", "red","white","blue"), values = scales::rescale((c(-9.5,-4,0,3.5))), limits= c(-9.5,3.5), name = "")+
  theme_minimal(base_size = 12)+
  theme(axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 10), panel.grid = element_blank())+
  labs(x="",y="")

dim(beta)<-355
data$beta<-beta
ggplot(data, aes(x=Tau, y=Commodity, fill=beta))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = c("red","white","blue"), values = scales::rescale((c(-2,0,2.7))), limits= c(-2,2.7), name = "")+
  theme_minimal(base_size = 12)+
  theme(axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 10), panel.grid = element_blank())+
  labs(x="",y="")

dim(lambda)<-355
data$lambda<-lambda
ggplot(data, aes(x=Tau, y=Commodity, fill=lambda))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = c("red2","red","white","blue","blue3"), values = scales::rescale((c(-6.1,-2,0,4,12.5))), limits= c(-6.1,12.5), name = "")+
  theme_minimal(base_size = 12)+
  theme(axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 10), panel.grid = element_blank())+
  labs(x="",y="")
##########################################################
coef01<-read.csv("coefrobust2-tau01.csv",header = TRUE)
coef025<-read.csv("coefrobust2-tau025.csv",header = TRUE)
coef05<-read.csv("coefrobust2-tau05.csv",header = TRUE)
coef075<-read.csv("coefrobust2-tau075.csv",header = TRUE)
coef09<-read.csv("coefrobust2-tau09.csv",header = TRUE)
pvalue01<-read.csv("pvaluerobust2-tau01.csv", header = TRUE)
pvalue025<-read.csv("pvaluerobust2-tau025.csv", header = TRUE)
pvalue05<-read.csv("pvaluerobust2-tau05.csv", header = TRUE)
pvalue075<-read.csv("pvaluerobust2-tau075.csv", header = TRUE)
pvalue09<-read.csv("pvaluerobust2-tau09.csv", header = TRUE)
coefalpha1<-cbind(coef01[,3],coef025[,3],coef05[,3],coef075[,3],coef09[,3])
coefbeta<-cbind(coef01[,4],coef025[,4],coef05[,4],coef075[,4],coef09[,4])
coeflambda<-cbind(coef01[,5],coef025[,5],coef05[,5],coef075[,5],coef09[,5])
pvaluealpha1<-cbind(pvalue01[,3],pvalue025[,3],pvalue05[,3],pvalue075[,3],pvalue09[,3])
pvaluebeta<-cbind(pvalue01[,4],pvalue025[,4],pvalue05[,4],pvalue075[,4],pvalue09[,4])
pvaluelambda<-cbind(pvalue01[,5],pvalue025[,5],pvalue05[,5],pvalue075[,5],pvalue09[,5])

alpha1<-matrix(NA,71,5)
for (i in 1:71) {
  for (j in 1:5) {
    if(pvaluealpha1[i,j]<0.1){
      alpha1[i,j]<-coefalpha1[i,j]
    }else{alpha1[i,j]<-0}
  }
}
beta<-matrix(NA,71,5)
for (i in 1:71) {
  for (j in 1:5) {
    if(pvaluebeta[i,j]<0.1){
      beta[i,j]<-coefbeta[i,j]
    }else{beta[i,j]<-0}
  }
}
lambda<-matrix(NA,71,5)
for (i in 1:71) {
  for (j in 1:5) {
    if(pvaluelambda[i,j]<0.1){
      lambda[i,j]<-coeflambda[i,j]
    }else{lambda[i,j]<-0}
  }
}


library(ggplot2)
library(dplyr)

commodities<-c("Corn, China", "Wheat, China", "Japonica Rice, China", "Pork, China", "Beef, China", 
               "Mutton, China", "Cotton, China", "Soybean, China", "Soybean Meal, China", 
               "Soybean Oil, China", "Rapeseed Meal, China", "Rapeseed Oil, China", "Palm Oil, China",
               "Pure Benzene, China", "Toluene, China", "O-Xylene, China", "Styrene, China", 
               "Butadiene, China", "P-Xylene, China", "PTA, China", "Epoxy Ethane, China", "Propylene, China", 
               "Polystyrene, China", "Natural Rubber, China", "PP, China", "Causticsoda, China", 
               "Potassic Fertilizers, China", "Xylene, China", "Glycol, China", "Acrylamide, China", 
               "Acrylic Acid, China", "Acrylonitrile, China", "Butyl Acrylate, China", "PET Chips, China", 
               "Methyl Methacrylate, China", "PVC, China", "Urea, China", "Methanol, China", 
               "Metallurgical Coke, China", "Hard Coking Coal, China", "Coking Coal, China", 
               "Hot-rolled Sheet/Coil, China", "Platinum, China", "Silver, China", "Palladium, China",
               "Gold, China", "Copper, China", "Nickel, China", "Aluminum, China", "Zinc, China",
               "Lead, China", "Tin, China", "Chrome Metal, China", "Silicon Metal, China", "Cobalt, China",
               "Ingot, China", "Electrolytic Manganese, China", "Nickel Sulfate, China", "Antimony, China",
               "Soybean Oil, US", "Diesel Oil, Singapore", "Gasoline, Japan", "Brent", "WTI", "Henry Hub Natural Gas",
               "Crude Oil, UAE Dubai", "Naphtha, Korea", "Gold, UK", "Silver, UK", "Platinum, UK", "Palladium, UK")
taus<-c("tau=0.1", "tau=0.25", "tau=0.5", "tau=0.75", "tau=0.9")
as.factor(taus)
data<-expand.grid(Commodity=commodities,Tau=taus)

dim(alpha1)<-355
data$alpha1<-alpha1
ggplot(data, aes(x=Tau, y=Commodity, fill=alpha1))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = c("red3", "red","white","blue"), values = scales::rescale((c(-7.3,-2,0,2.2))), limits= c(-7.3,2.2), name = "")+
  theme_minimal(base_size = 12)+
  theme(axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 10), panel.grid = element_blank())+
  labs(x="",y="")

dim(beta)<-355
data$beta<-beta
ggplot(data, aes(x=Tau, y=Commodity, fill=beta))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = c("red","white","blue","blue3"), values = scales::rescale((c(-3,0,3,10))), limits= c(-3,10), name = "")+
  theme_minimal(base_size = 12)+
  theme(axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 10), panel.grid = element_blank())+
  labs(x="",y="")

dim(lambda)<-355
data$lambda<-lambda
ggplot(data, aes(x=Tau, y=Commodity, fill=lambda))+
  geom_tile(color="white")+
  scale_fill_gradientn(colors = c("red2","red","white","blue","blue2"), values = scales::rescale((c(-4.5,-2,0,3,5))), limits= c(-4.5,5), name = "")+
  theme_minimal(base_size = 12)+
  theme(axis.text.y=element_text(size = 7), axis.text.x=element_text(size = 10), panel.grid = element_blank())+
  labs(x="",y="")
