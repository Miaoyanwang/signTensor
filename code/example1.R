source("signT.R")
library(ggplot2)
set.seed(1)
d=30;r=3;
data_full=NULL
for(nsim in 1:10){
a=matrix(rnorm(d*r),ncol=r)
T=tensorize(a,a,a)
c=c(0.1,0.5,1,3,5,10,15,20,25,30,35,40,45,50,75,100,125,150,175,200)
rank2=accuracy=NULL
for(i in 1:20){
T_trans=1/(1+exp(-(c[i]*T)))
T_trans=as.tensor(T_trans/sqrt(sum(T_trans^2)))
tensor_rank=appx_rank(T_trans,thresh=90,step=5)
rank2=c(rank2,max(tensor_rank[,1]))
accuracy=c(accuracy,tensor_rank[dim(tensor_rank)[1],2])
}
data=data.frame(c,rank2)
colnames(data)=c("c","rank")
data_full=rbind(data_full,data)
}
data_full=cbind(data_full,rep(1:10,rep(20,10)))
colnames(data_full)=c("c","rank","sim")
save(data_full,file="example1.RData")

ag <- aggregate(. ~ c, data_full, function(x) c(mean = mean(x), sd = sd(x)))

plotc=c(1:8,10,12,14:20)
ag=ag[plotc,]
pdf("example1.pdf",width=3.5,height=3.1)
ggplot(ag,aes(c,rank[,1]))+geom_point(aes(c,rank[,1]))+geom_line(aes(c,rank[,1]))+geom_errorbar(aes(ymin=rank[,1]-rank[,2]/sqrt(10), ymax=rank[,1]+rank[,2]/sqrt(10)), width=.2,position=position_dodge(.9))+labs(x = "transformation level (c)")+labs(y = "numerical rank")
dev.off()
