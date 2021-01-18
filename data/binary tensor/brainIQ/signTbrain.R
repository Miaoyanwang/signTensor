args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}
source("signT.R")
load(file = "brain_binary_IQ.RData")




rk = rep(c(3,6,9,12,15),each = 5)


tnsr = A
nfoldlist = c(5,3,2,3)
indexset = expand.grid(rk,nfoldlist)
r = indexset[BATCH,1]
nfold = indexset[BATCH,2]


#kfoldindex = c("25missing","33missing","50missing","67missing")[(BATCH-1)%/%25+1]
set.seed(BATCH)
#cindex = split(sample(1:length(tnsr),length(tnsr)),as.factor(1:5))
l1 = split(sample(which(tnsr==1),length(which(tnsr==1))),as.factor(1:nfold))
l2 = split(sample(which(tnsr==0),length(which(tnsr==0))),as.factor(1:nfold))

cindex = list()
for (k in 1:nfold) {
  cindex[[k]] = c(l1[[k]],l2[[k]])
}


cvr = vector(length = nfold)
cvrr = vector(length = nfold)
cvr2 = vector(length = nfold)
cvrr2 = vector(length = nfold)

# Est = list()
# Est2 = list()
for(i in 1:nfold){
  if (BATCH<76) {
    testindex =  cindex[[i]]
  }else{
    testindex = setdiff(1:length(tnsr),cindex[[i]])
  }
  #testindex = cindex[[i]] # 33 percent available
  train_T = tnsr
  train_T[testindex] = NA
  res=SignT(train_T,r,Lmin=0,Lmax=1,H=20,option=2) ## recommend option = 2 or 3.
  res2 = fit_continuous(train_T,r)
  
  # Est[[i]] = res$est
  # Est2[[i]] = res2$est
  cvr[i] = mean(abs(res$est[testindex]-tnsr[testindex]))
  cvrr[i] = mean(abs(round(res$est[testindex])-tnsr[testindex]))
  cvr2[i] = mean(abs(res2$est[testindex]-tnsr[testindex]))
  cvrr2[i] = mean(abs(round(res2$est[testindex])-tnsr[testindex]))
  
}

save(cvrr2,cvrr,cvr2,cvr,file = paste("signT_braincv2_",BATCH,".RData",sep = ""))


###############################################################################################
##### MAE cross validation result when r = 3,6,9,12,15 with different missing rate(20,33,50,67)
###############################################################################################
mcr_cv= matrix(nrow = 20,ncol = 6)
mcr_cv[,6] = rep(c(20,33,50,67),each = 5)
mcr_cv2 = matrix(nrow = 20,ncol = 6)
mcr_cv2[,6] = rep(c(20,33,50,67),each = 5)

for (BATCH in 1:100) {
  load(file = paste("CV/signT_braincv2_",BATCH,".RData",sep = ""))

  row_ind = BATCH%%5+1+5*((BATCH-1)%/%25)
  col_ind = ((BATCH-1)%%25)%/%5+1
  mcr_cv[row_ind,col_ind] = mean(cvr)
  mcr_cv2[row_ind,col_ind] = mean(cvr2)

}



meancv = matrix(nrow = 4,ncol = 5)
sdcv = matrix(nrow  = 4,ncol = 5)

meancv2 = matrix(nrow = 4,ncol = 5)
sdcv2 = matrix(nrow  = 4,ncol = 5)

for(i in 1:5){
  for(j in 1:4){
    meancv[j,i] = mean(mcr_cv[,i][(1+5*(j-1)):(5*j)])
    sdcv[j,i] = sd(mcr_cv[,i][(1+5*(j-1)):(5*j)])
    meancv2[j,i] = mean(mcr_cv2[,i][(1+5*(j-1)):(5*j)])
    sdcv2[j,i] = sd(mcr_cv2[,i][(1+5*(j-1)):(5*j)])
  }

}

###########################################################################################
##### Rank 1:60, 4 different missing rates (20,33,50,67)
##########################################################################################

mcr_cv = NULL
mcr_cv2 = NULL
auc_cv = NULL
auc_cv2 = NULL
fnr_cv = NULL
fnr_cv2 = NULL
fpr_cv = NULL
fpr_cv2 = NULL
for (i in 1:60) {
  r = i
  load(file = paste("CV/signT_braincv2_",r,".RData",sep = ""))
  a = 0; a2 = 0
  n = 0; n2 =0
  p = 0; p2 = 0
  for(j in 1:3){
    a = a+auc(tnsr[testindex],round(Est[[j]][testindex]))
    a2 = a2+auc(tnsr[testindex],round(Est2[[j]][testindex]))
    n = n +  length(which(tnsr[testindex]-round(Est[[j]][testindex])==1)) /length(which(tnsr[testindex]==0))
    n2 = n2 + length(which(tnsr[testindex]-round(Est2[[j]][testindex])==1)) /length(which(tnsr[testindex]==0))
    p = p +  length(which(tnsr[testindex]-round(Est[[j]][testindex])==-1)) /length(which(tnsr[testindex]==1))
    p2 = p2 + length(which(tnsr[testindex]-round(Est2[[j]][testindex])==-1)) /length(which(tnsr[testindex]==1))
  }

  auc_cv = c(auc_cv,a/5)
  auc_cv2 = c(auc_cv2,a2/5)
  fnr_cv = c(fnr_cv,n/5)
  fnr_cv2 = c(fnr_cv2,n2/5)
  fpr_cv = c(fpr_cv,p/5)
  fpr_cv2 = c(fpr_cv2,p2/5)

  mcr_cv = c(mcr_cv,mean(cvr))
  mcr_cv2 = c(mcr_cv2,mean(cvr2))
}

par(mfrow = c(2,2))
plot(1:60,mcr_cv,xlab = "rank",main = "MCR result",pch =20,type = "l",col ="red")
points(1:60,mcr_cv2,pch = 20,type ="l")
legend(40, 0.12, legend=c("NonParaT", "CPT"),
       col=c( "red","black"),lty = c(1,1))

plot(1:60,auc_cv,xlab = "rank",main = "AUC result",pch =20,type = "l",col ="red")
points(1:60,auc_cv2,pch = 20,type = "l")
legend(40, 0.87, legend=c("NonParaT", "CPT"),
       col=c("red","black"),lty = c(1,1))

plot(1:60,fnr_cv,xlab = "rank",main = "FNR result",pch =20,type = "l",col = "red")
points(1:60,fnr_cv2,pch = 20,type= "l")
legend(40, 0.11, legend=c("NonParaT", "CPT"),
       col=c("red","black"),lty = c(1,1))

plot(1:60,fpr_cv,xlab = "rank",main = "FPR result",pch =20,type = "l",col ="red")
points(1:60,fpr_cv2,pch = 20,type="l")
legend(40, 0.12, legend=c("NonParaT", "CPT"),
       col=c("black", "red"),lty  = c(1,1))


###########################################################################################
#### Brain visualization
#### output with rank =10, H = 20

load(file = "CV/signT_brain_2.RData")

Est = res$est

numH = length(FSIQ[FSIQ>quantile(FSIQ,0.75)])
Hbrain = Est[,,which(FSIQ>quantile(FSIQ,0.75))]

AvgHbrain = 0
for(i in 1:numH){
  AvgHbrain = AvgHbrain + Hbrain[,,i]
}
AvgHbrain = (AvgHbrain+t(AvgHbrain))/(2*numH)

image(AvgHbrain)


numL = length(FSIQ[FSIQ<quantile(FSIQ,0.25)])
Lbrain = Est[,,which(FSIQ<quantile(FSIQ,0.25))]
AvgLbrain = 0
for(i in 1:numL){
  AvgLbrain = AvgLbrain + Lbrain[,,i]
}
AvgLbrain = (AvgLbrain+t(AvgLbrain))/(2*numL)
image(AvgLbrain)
max(AvgHbrain-AvgLbrain)

diffbrain = AvgHbrain-AvgLbrain

image(diffbrain)
index = which(diffbrain>0.11,arr.ind = T)

diffedge = matrix(0,nrow = 68,ncol = 68)
for(i in 1:20){
  diffedge[index[i,][1],index[i,][2]]=1
}

write(diffedge,file = "diff.edge")
node = read.table("/Users/chanwoolee/Desktop/Study/Research/BrainNetViewer_20171031/Data/miaoyan/Desikan_miaoyan.node")
node2 = read.table("/Users/chanwoolee/Desktop/Study/Research/BrainNetViewer_20171031/Data/miaoyan/ROI_68_lobe_degree.node")
node[,6]


head(node)
head(node2)
