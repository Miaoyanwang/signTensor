### min hypergraphon, type=9
source("signT.R")
set.seed(1)
error=error_con=error_m=NULL
d_list=seq(from=15,to=60,length=10)

args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
    BATCH <- as.numeric(args[1])
} else {
    stop()
}

truer=2
d=d_list[BATCH]
##### simulate graphon models
a=seq(from=0,to=1,length=d)
signal=graphon_to_tensor(a,a,a,type=9)
for(sim in 1:30){
Y=signal+array(rnorm(length(signal),0,0.15),dim=dim(signal)) ###
Lmin=min(Y,na.rm=T)
Lmax=max(Y,na.rm=T)
set.seed(1)
res=SignT(Y,truer,Lmin=Lmin,Lmax=Lmax,H=10+2*(d-15)/10,option=2)
error=c(error,mean(abs(res$est-signal)))
### matrix method
res_m=SignUnfold(Y,truer,Lmin=Lmin,Lmax=Lmax,H=10+2*(d-15)/10,option=2,signal)
error_m=c(error_m,res_m$error)
##### continuous method
res2=fit_continuous(Y,truer)
error_con=c(error_con,mean(abs(res2$est-signal)))
print(paste("dimension",d,"-simulation",sim," is done",sep = ""))
}

save(error,error_con,error_m,file=paste("Figure4-dimension",d,".RData",sep=""))


####################### load output from server#######################
d_list=seq(from=15,to=60,length=10)
con=ours=NULL
load("Figure4-dimension15.RData")
ours=error; con=error_con; matrix=error_m
for(i in d_list[2:10]){
    load(sprintf("Figure4-dimension%d.RData",i))
    ours=rbind(ours,error)
    con=rbind(con,error_con)
    matrix=rbind(matrix,error_m)
}
error=ours;error_con=con;error_m=matrix
save(error,error_con,error_m,file="Model4.RData")
#######################


##### plot##########
d_list=seq(from=15,to=60,length=10)
load("Model4.RData")
library(ggplot2)
error_cont=cbind(d_list,apply(error_con,1,mean),apply(error_con,1,sd))
error_ours=cbind(d_list,apply(error,1,mean),apply(error,1,sd))
error_m=cbind(d_list,apply(error_m,1,mean),apply(error_m,1,sd))
data=data.frame(rbind(error_cont,error_ours,error_m))
data=cbind(data,c(rep("cont",10),rep("ours",10),rep("matrix",10)))
colnames(data)=c("dim","mean","se","method")

pdf("Model4.pdf",width=4.5,height=3)
figure=ggplot(data,aes(dim,mean))+geom_point(aes(dim,mean,col=method,shape=method),size=2)+geom_line(aes(dim,mean,col=method))+geom_errorbar(aes(ymin=mean-se/sqrt(30), ymax=mean+se/sqrt(30),col=method), width=1,position=position_dodge(0))+labs(x = "dimension")+labs(y = "MAE")

figure
dev.off()

