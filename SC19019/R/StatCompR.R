#knn estimate function for standard normal distribution,find the best k
library(grf)
library(MASS)
options(digits=6)
m=50
s=3
p=11
nn=0
generate_fixed_point<-function(j)
{
  x<-rep(1,j)
  x<-x-1
  for(i in 3:j)
  {
    q=(i-1)%/%2
    x[q]=0.5
  }
  return(x)
}
generate_sample<-function(nn,m,p)
{
  set.seed(1234+nn)
  mean<-rep(1,p)
  sigma<-diag(mean)
  mean<-mean-1
  mydata<-mvrnorm(m,mean,sigma)
  w<-rbinom(m,1,0.5)
  z=cbind(mydata,w)
  return(z)
}

generate_setting1_2<-function(mydata,m,p)
{
  x<-mydata[,1:(p-1)]
  w<-mydata[,(p+1)]
  if(p==4)
    y=(x[,1]-1)^2*w+(x[,2]+1)^3*w-3*x[,3]*w+mydata[,4]
  if(p==11)
    y=(x[,3]-1)^2*w+(x[,5]+1)^3*w-3*x[,7]*w+mydata[,11]
  z=cbind(mydata,y)
  return(z)
}

generate_setting3_11<-function(mydata,m,j)
{
  y<-rep(1,m)
  y=y-1
  w<-mydata[,(j+1)]
  for(i in 1:(j-1))
  {
    y=y+((mydata[,i]^3-2*mydata[,i]^2+2*mydata[,i])^2)
  }
  y=log(y)*w+mydata[,j]
  z=cbind(mydata,y)
  return(z)
}
estimate_setting_DNN<-function(zz,mm,ss,x,p)
{
  d<-matrix(nrow=mm,ncol=2)
  for(i in 1:mm)
  {
    d[i,1]=t(zz[i,1:(p-1)]-x)%*%(zz[i,1:(p-1)]-x)
  }
  d[,2]=zz[,(p+2)]
  pp=d[order(d[,1]),]	
  DNN=0
  for(i in 1:(mm-ss+1))
  {
    k=choose((mm-i),(ss-1))
    h=choose(mm,ss)
    DNN=DNN+k/h*pp[i,2]
  }
  
  return(DNN)
}

twoscaled_DNN<-function(zz,p,s,mm,xx)
{
  z2<-data.frame(zz)
  data1<-subset(x=z2,z2$w==1)
  data2<-subset(x=z2,z2$w==0)
  new1<-as.matrix(data1)
  new2<-as.matrix(data2)
  dd1<-dim(new1)
  dd2<-dim(new2)
  mm1=dd1[1]
  mm2=dd2[1]
  dnn1_s=estimate_setting_DNN(new1,mm1,s,xx,p)
  ss=2*s
  dnn1_2s=estimate_setting_DNN(new1,mm1,ss,xx,p)
  kk1=1/(1-2^(2/(p-1)))
  kk2=-2^(2/(p-1))*kk1
  twoscales_dnn1=kk1*dnn1_s+kk2*dnn1_2s
  dnn2_s=estimate_setting_DNN(new2,mm2,s,xx,p)
  dnn2_2s=estimate_setting_DNN(new2,mm2,ss,xx,p)
  twoscales_dnn2=kk1*dnn2_s+kk2*dnn2_2s
  
  effect=twoscales_dnn1-twoscales_dnn2
  return(effect)
}

generate_CF_sample<-function(mydata,m,p,x.test)
{
  x<-mydata[,1:(p-1)]
  w<-mydata[,(p+1)]
  y<-mydata[,(p+2)]
  c.forest=causal_forest(x,y,w)
  c.pred=predict(c.forest,t(x.test))
  c.pre=c.pred[1,1]
  return(c.pre)
}
result_setting1_2<-function(m,s,p)
{
  estimate_result<-matrix(0,nrow=2,ncol=6)
  two_DNN<-vector("numeric",m)
  cf_mm<-vector("numeric",m)
  if(p==4)
    x.test<-c(0.5,-0.5,0.5)
  if(p==11)
    x.test<-c(0,0,0.5,0,-0.5,0,0.5,0,0,0)
  for(i in 1:m)
  {
    mydata=generate_sample(i,m,p)
    zz=generate_setting1_2(mydata,m,p)
    two_DNN[i]=twoscaled_DNN(zz,p,s,m,x.test)
    cf_mm[i]=generate_CF_sample(zz,m,p,x.test)
  } 
  if(p==4)
  {    
    y=(x.test[1]-1)^2+(x.test[2]+1)^3-3*x.test[3]
  } 
  if(p==11)
  {    
    y=(x.test[3]-1)^2+(x.test[5]+1)^3-3*x.test[7]
  } 
  
  estimate_result[1,1]=mean(two_DNN)-y
  estimate_result[2,1]=mean(cf_mm)-y
  estimate_result[1,2]=t(two_DNN-y)%*%(two_DNN-y)/m
  estimate_result[2,2]=t(cf_mm-y)%*%(cf_mm-y)/m
  estimate_result[1,3]=var(two_DNN)
  estimate_result[2,3]=var(cf_mm)
  
  two_DNN2<-vector("numeric",m)
  cf_mm2<-vector("numeric",m)
  x2.test<-runif((p-1), min = 0, max = 1)
  for(i in 2:(m+1))
  {
    mydata=generate_sample(i,m,p)
    zz=generate_setting1_2(mydata,m,p)
    two_DNN2[i]=twoscaled_DNN(zz,p,s,m,x2.test)
    cf_mm2[i]=generate_CF_sample(zz,m,p,x2.test)
  } 
  
  if(p==4)
    y2=(x2.test[1]-1)^2+(x2.test[2]+1)^3-3*x2.test[3]
  if(p==11)
    y2=(x2.test[3]-1)^2+(x2.test[5]+1)^3-3*x2.test[7]
  
  estimate_result[1,5]=mean(two_DNN2)-y2
  estimate_result[2,5]=mean(cf_mm2)-y2
  estimate_result[1,6]=t(two_DNN2-y2)%*%(two_DNN2-y2)/m
  estimate_result[2,6]=t(cf_mm2-y2)%*%(cf_mm2-y2)/m
  return(estimate_result)
}
generate_response3_11<-function(x,j)
{
  
  y=0
  for(i in 1:(j-1))
  {
    y=y+((x[i]^3-2*x[i]^2+2*x[i])^2)
  }
  y=log(y)
  return(y)
  
}
result_setting3_11<-function(m,s,p)
{
  estimate_result<-matrix(0,nrow=2,ncol=6)
  two_DNN<-vector("numeric",m)
  cf_mm<-vector("numeric",m)
  x.test=generate_fixed_point((p-1))
  for(i in 2:(m+1))
  {
    mydata=generate_sample(i,m,p)
    #head(mydata,n=5)
    zz=generate_setting3_11(mydata,m,p)
    #head(zz,n=5)
    two_DNN[i]=twoscaled_DNN(zz,p,s,m,x.test)
    cf_mm[i]=generate_CF_sample(zz,m,p,x.test)
  } 
  y=generate_response3_11(x.test,(p-1))
  estimate_result[1,1]=mean(two_DNN)-y
  # cfmm=as.vector(cf_mm)
  estimate_result[2,1]=mean(cf_mm)-y
  estimate_result[1,2]=t(two_DNN-y)%*%(two_DNN-y)/m
  estimate_result[2,2]=t(cf_mm-y)%*%(cf_mm-y)/m
  estimate_result[1,3]=var(two_DNN)
  estimate_result[2,3]=var(cf_mm)
  x2.test<-runif((p-1), min = 0, max = 1)
  two_DNN2<-vector("numeric",m)
  cf_mm2<-vector("numeric",m)
  for(i in 1:m)
  {
    mydata=generate_sample(i,m,p)
    zz=generate_setting3_11(mydata,m,p)
    two_DNN2[i]=twoscaled_DNN(zz,p,s,m,x2.test)
    cf_mm2[i]=generate_CF_sample(zz,m,p,x2.test)
  } 
  
  
  y2=generate_response3_11(x2.test,(p-1))
  estimate_result[1,5]=mean(two_DNN2)-y2
  estimate_result[2,5]=mean(cf_mm2)-y2
  estimate_result[1,6]=t(two_DNN2-y2)%*%(two_DNN2-y2)/m
  estimate_result[2,6]=t(cf_mm2-y2)%*%(cf_mm2-y2)/m
  return(estimate_result)
}
summary_result<-function(m,s)
{
  re1=result_setting1_2(m,s,4)
  re2=result_setting1_2(m,s,11)
  re=rbind(re1,re2)
  re3=result_setting3_11(m,s,11)
  re=rbind(re,re3)
  re4=result_setting3_11(m,s,16)
  re=rbind(re,re4)
  re5=result_setting3_11(m,s,21)
  re=rbind(re,re5) 
  re6=result_setting3_11(m,s,26)
  re=rbind(re,re6) 
  re7=result_setting3_11(m,s,31)
  re=rbind(re,re7) 
  re8=result_setting3_11(m,s,36)
  re=rbind(re,re8) 
  re9=result_setting3_11(m,s,41)
  re=rbind(re,re9) 
  re10=result_setting3_11(m,s,45)
  re=rbind(re,re10) 
  re11=result_setting3_11(m,s,51)
  re=rbind(re,re11) 
  return(re)
}


