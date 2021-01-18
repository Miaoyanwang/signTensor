source("signT.R")
library(R.matlab)
## test
set.seed(1)
d=5

##### simulate graphon models
a=seq(from=0,to=1,length=d)
b=seq(from=0,to=1,length=d)
c=seq(from=0,to=1,length=d)
signal=graphon_to_tensor(a,b,c,type=10)
Y=signal+array(rnorm(length(signal),0,0*max(abs(signal))),dim=dim(signal))
truer=2
hist(Y)
missing=array(rbinom(length(signal),1,0),dim=dim(signal))
Y[missing==1]=NA
Lmin=min(Y,na.rm=T)
Lmax=max(Y,na.rm=T)

set.seed(1)
res=SignT(Y,truer,Lmin=Lmin,Lmax=Lmax,H=10,option=2) ## recommend option = 2 or 3
plot(res$est,signal)
abline(0,1)
plot(res$est[missing==0],signal[missing==0])
plot(res$est[missing==1],signal[missing==1])
mean(abs(res$est[missing==0]-signal[missing==0]))
abline(0,1)


### continuous
res2=fit_continuous(Y,truer)
plot(res2$est,signal)
plot(res2$est[missing==0],signal[missing==0])
plot(res2$est[missing==1],signal[missing==1])
mean(abs(res2$est[missing==0]-signal[missing==0]))
abline(0,1)

data=readMat("../data/dnations.mat")
#data=readMat("../data/alyawarradata.mat")
tensor=data$R-mean(data$R,na.rm=T)
set.seed(1)

hold=sample(sum(is.na(tensor)==F),0.2*sum(is.na(tensor)==F),replace=TRUE) ## test
training=tensor
training[hold]=NA

truer=3
res=SignT(training,truer,min(training,na.rm=T),max(training,na.rm=T),H=20,option=2)
res2=fit_continuous(training,truer)

mean(abs(res$est[hold]-tensor[hold]),na.rm=T) ## ours
mean(abs(res2$est[hold]-tensor[hold]),na.rm=T)

library("pROC")
auc(tensor[hold],res$est[hold]) ## ours
auc(tensor[hold],res2$est[hold])

