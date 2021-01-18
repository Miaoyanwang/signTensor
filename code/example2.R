source("signT.R")
library(ggplot2)
set.seed(1)
d=30;
Z=T=array(0,dim=c(d,d,d))
for(i in 1:d){
    for(j in 1:d){
        for(k in 1:d){
        Z[i,j,k]=max(c(i,j,k))/d
        T[i,j,k]=log(1+Z[i,j,k])
    }
}
}
res=cp(as.tensor(Z),30)
res2=cp(as.tensor(T),30)

data=cbind(1:30,sort(res$lambda,decreasing=TRUE))
data=rbind(data,cbind(1:30,sort(res2$lambda,decreasing=TRUE)))
colnames(data)=c("component","singular value")
tensor=as.factor(c(rep("Z",30),rep("T",30)))
data[,2]=log(data[,2])
data=data.frame(data,tensor)

pdf("example2.pdf",width=4,height=3)
ggplot(data,aes(component,singular.value))+geom_point(aes(component,singular.value,shape=tensor,col=tensor))+ylim(6,10.1)+labs(x = "index")+labs(y = "singular value on log scale")
dev.off()
