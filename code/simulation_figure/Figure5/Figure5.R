## Model 5. stochastic tensor block
source("signT.R")
set.seed(1)
error=error_con=error_m=NULL
d=40
rho_list=seq(from=0.3,to=1,length=8)

args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
    BATCH <- as.numeric(args[1])
} else {
    stop()
}

truer=3
rho=rho_list[BATCH]
obs=rbinom(d^3,1,rho)
##### simulate graphon models
a=seq(from=0,to=1,length=d)
signal=graphon_to_tensor(a,a,a,type=5)
for(sim in 1:30){
Y=signal+array(runif(length(signal),-0.3,0.3),dim=dim(signal))
Y[obs==0]=NA
Lmin=min(Y,na.rm=T)
Lmax=max(Y,na.rm=T)
##### our method
res=SignT(Y,truer,Lmin=Lmin,Lmax=Lmax,H=10+2*(d-15)/10,option=2)
error=c(error,mean(abs(res$est-signal)))
#### matrix method
res_m=SignUnfold(Y,truer,Lmin=Lmin,Lmax=Lmax,H=10+2*(d-15)/10,option=2,signal)
error_m=c(error_m,res_m$error)
##### continuous method
res2=fit_continuous(Y,truer)
error_con=c(error_con,mean(abs(res2$est-signal)))
print(paste("dimension",d,"-simulation",sim," is done",sep = ""))
}

save(error,error_con,error_m,file=paste("Figure5-dimension",rho,".RData",sep=""))


####################### load output from server#######################
rho_list=seq(from=0.3,to=1,length=8)
load("Figure5-dimension0.3.RData")
ours=error; con=error_con; matrix=error_m
for(i in 2:7){
    load(sprintf("Figure5-dimension%.1f.RData",rho_list[i]))
    ours=rbind(ours,error)
    con=rbind(con,error_con)
    matrix=rbind(matrix,error_m)
}
load("Figure5-dimension1.RData")
ours=rbind(ours,error); con=rbind(con,error_con); matrix=rbind(matrix,error_m)
error=ours;error_con=con;error_m=matrix
save(error,error_con,error_m,file="Model5.RData")
#######################

##### plot
rho_list=seq(from=0.3,to=1,length=8)
load("Model5.RData")
library(ggplot2)
error_cont=cbind(rho_list,apply(error_con,1,mean),apply(error_con,1,sd))
error_ours=cbind(rho_list,apply(error,1,mean),apply(error,1,sd))
error_m=cbind(rho_list,apply(error_m,1,mean),apply(error_m,1,sd))
data=data.frame(rbind(error_cont,error_ours,error_m))
data=cbind(data,c(rep("cont",8),rep("ours",8),rep("matrix",8)))
colnames(data)=c("dim","mean","se","method")


pdf("Model5.pdf",width=4.5,height=3)
figure=ggplot(data,aes(x=dim,y=mean))+geom_point(aes(x=dim,y=mean,col=method,shape=method),size=2)+geom_line(aes(x=dim,y=mean,col=method))+geom_errorbar(aes(ymin=mean-se/sqrt(30), ymax=mean+se/sqrt(30),col=method), width=.01,position=position_dodge(0))+labs(x = "Observation fraction")+labs(y = "MAE")
#+coord_cartesian(ylim = c(0.02, 0.16))

figure
dev.off()

