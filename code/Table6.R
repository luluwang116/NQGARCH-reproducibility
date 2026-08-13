return<-read.csv("final return 2016.csv",header = TRUE)
return<-return[,-1]
return0<-ts(return)
dim(return0)
return0<-t(return0)
library(moments)
library(tseries)
stat<-matrix(0,71,8)
for (i in 1:71) {
  stat[i,1]<-mean(return0[i,])
  stat[i,2]<-sd(return0[i,])
  stat[i,3]<-min(return0[i,])
  stat[i,4]<-max(return0[i,])
  stat[i,5]<-skewness(return0[i,])
  stat[i,6]<-kurtosis(return0[i,])
  adf<-adf.test(return0[i,])
  stat[i,7]<-adf$p.value
  jb<-jarque.bera.test(return0[i,])
  stat[i,8]<-jb$p.value
}
round(stat,4)