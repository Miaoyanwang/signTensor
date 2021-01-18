library(rTensor)
library(quadprog)
library(Matrix)
min.thresh=10^(-3)
max.thresh=10^10

binaryloss=function(Ybar,W,Yfit){
    return(mean(W*abs(Ybar-sign(Yfit)),na.rm=TRUE))
}
SignUnfold=function(Y,truer,Lmin=Lmin,Lmax=Lmax,H=10,option,signal){
    res=res1=SignT(array(unfold(as.tensor(Y),2:3,1)@data,dim=c(d^2,d,1)),truer,H=H,Lmin=Lmin,Lmax=Lmax,option)
    error=error_1=mean(abs(res1$est-array(unfold(as.tensor(signal),2:3,1)@data,dim=c(d^2,d,1))))
     return(list("est"=res,"error"=error))
     
    res2=SignT(array(unfold(as.tensor(Y),c(1,3),2)@data,dim=c(d^2,d,1)),truer,H=H,Lmin=Lmin,Lmax=Lmax,option)
    error_2=mean(abs(res2$est-array(unfold(as.tensor(signal),c(1,3),2)@data,dim=c(d^2,d,1))))
    if(error_2<error){res=res2; error=error_2}
   
   res3=SignT(array(unfold(as.tensor(Y),1:2,3)@data,dim=c(d^2,d,1)),truer,H=H,Lmin=Lmin,Lmax=Lmax,option)
error_3=mean(abs(res3$est-array(unfold(as.tensor(signal),1:2,3)@data,dim=c(d^2,d,1))))
     if(error_3<error){res=res3; error=error_3}
   
    return(list("est"=res,"error"=(error_1+error_2+error_3)/3))
}
#################### main function for nonparametric tensor completion  ####################
SignT=function(Y,truer,H=5,Lmin,Lmax,rho=0.1,lambda=10^(-3),option=2){
    B_fitted=result=list()
    pi_seq=seq(from=Lmin,to=Lmax,length=2*H+1)
    
    for(h in 2:(2*H)){
        pi=pi_seq[h]
        if(option==1){
            res=ADMM(sign(Y-pi),abs(Y-pi),r=truer,rho=rho,lambda=lambda)}
        else if(option==2){
            res=Alt(sign(Y-pi),abs(Y-pi),r=truer,type="logistic",start="linear")## recommend
            
        }else if(option==3){
            res=Alt(sign(Y-pi),abs(Y-pi),r=truer,type="hinge",start="linear")## recommend
        }
        result[[h]]=res
        B_fitted[[h]]=res$fitted
        print(paste("----",h,"-th level finished --- ",sep=""))
    }
    B_fitted=array(unlist(B_fitted),dim=c(dim(Y),2*H-1));
    res=list();
    res$result=result;
    res$fitted=B_fitted
    res$est=1/2*(apply(sign(B_fitted),1:3,sum)/(2*H)+1)*(Lmax-Lmin)+Lmin
    return(res)
}

### Alternating optimization for classification
Alt=function(Ybar,W,r,type=c("logistic","hinge"),start="random"){
    result=list()
    d=dim(Ybar)
    if(start=="linear"){
    sink("NULL")
    ini=fit_continuous(Ybar,r)
    
    #ini=fit_continuous(tensorize(a,b,c)-quantile(tensorize(a,b,c),mean((Ybar*W<0))),r);
    #diag(scale)=ini$lambda
    #A1=cbind(a,rep(1,d[1]))
    #A2=cbind(b,rep(1,d[2]))
    #A3=cbind(c,-quantile(tensorize(a,b,c),mean((Ybar*W<0)))*rep(1,d[3]))
    
    sink()
    A1 = ini$U[[1]];
    A2 = ini$U[[2]];
    scale=matrix(0,nrow=r,ncol=r)
    diag(scale)=ini$lambda
    A3 = ini$U[[3]]%*%scale;
    }else{
    A1 = cbind(randortho(d[1])[,1:r]);
    A2 = cbind(randortho(d[2])[,1:r]);
    A3 = cbind(randortho(d[3])[,1:r]);
    }
    obj=cost(A1,A2,A3,Ybar,W,type);
    binary_obj=binaryloss(Ybar,W,tensorize(A1,A2,A3))
    
    error=1;iter=1;
 
 while((error>0.01)&(binary_obj[iter]>0.01)&(iter<20)){
     
 
 #tic()
 optimization=optim(c(A3),function(x)cost(A1,A2,matrix(x,ncol=r),Ybar,W,type),function(x)gradient(A1,A2,matrix(x,ncol=r),3,Ybar,W,type),method="BFGS")
 A3=matrix(optimization$par,ncol=r)
 optimization=optim(c(A2),function(x)cost(A1,matrix(x,ncol=r),A3,Ybar,W,type),function(x)gradient(A1,matrix(x,ncol=r),A3,2,Ybar,W,type),method="BFGS")
 A2=matrix(optimization$par,ncol=r)
 optimization=optim(c(A1),function(x)cost(matrix(x,ncol=r),A2,A3,Ybar,W,type),function(x)gradient(matrix(x,ncol=r),A2,A3,1,Ybar,W,type),method="BFGS")
 A1=matrix(optimization$par,ncol=r)
 #toc()

 #tic()
 #l=lapply(1:nrow(A3),function(i){optim(c(A3[i,]),function(x)cost(A1,A2,matrix(x,ncol=r),array(Ybar[,,i],dim=c(d[1],d[2],1)),array(W[,,i],dim=c(d[1],d[2],1)),type),function(x)gradient(A1,A2,matrix(x,ncol=r),3,array(Ybar[,,i],dim=c(d[1],d[2],1)),array(W[,,i],dim=c(d[1],d[2],1)),type),method="BFGS")$par}) ## perhaps faster than the other approach?
 #A3=matrix(unlist(l),nrow=nrow(A3),byrow=T)
 #l=lapply(1:nrow(A2),function(i){optim(c(A2[i,]),function(x)cost(A1,matrix(x,ncol=r),A3,array(Ybar[,i,],dim=c(d[1],1,d[3])),array(W[,i,],dim=c(d[1],1,d[3])),type),function(x)gradient(A1,matrix(x,ncol=r),A3,2,array(Ybar[,i,],dim=c(d[1],1,d[3])),array(W[,i,],dim=c(d[1],1,d[3])),type),method="BFGS")$par}) ## perhaps faster than the other approach?
 #A2=matrix(unlist(l),nrow=nrow(A2),byrow=T)
 #l=lapply(1:nrow(A1),function(i){optim(c(A1[i,]),function(x)cost(matrix(x,ncol=r),A2,A3,array(Ybar[i,,],dim=c(1,d[2],d[3])),array(W[i,,],dim=c(1,d[2],d[3])),type),function(x)gradient(matrix(x,ncol=r),A2,A3,1,array(Ybar[i,,],dim=c(1,d[2],d[3])),array(W[i,,],dim=c(1,d[2],d[3])),type),method="BFGS")$par}) ### perhaps faster than the other approach?
 #A1=matrix(unlist(l),nrow=nrow(A1),byrow=T)
 #toc()
        obj=c(obj,cost(A1,A2,A3,Ybar,W,type))
        binary_obj=c(binary_obj,binaryloss(Ybar,W,tensorize(A1,A2,A3)))
        iter=iter+1
error=(obj[iter-1]-obj[iter])

 }
 result$binary_obj=binary_obj;
 result$obj=obj;
 result$iter=iter;
 result$error=error;
 result$fitted=tensorize(A1,A2,A3); ## exact low-rank
 return(result)
}

gradient=function(A1,A2,A3,mode,Ybar,W,type=c("logistic","hinge")){
    d=dim(Ybar)
    margin=Ybar*tensorize(A1,A2,A3)
    R=dim(A3)[2]
    
    if(type=="logistic"){
        tem=-W*Ybar*exp(-margin)/(1+exp(-margin))
    }else if(type=="hinge"){
        tem=-W*Ybar*(margin<1)
    }
    tem[is.na(tem)]=0
    
    if(mode==3){
    Grad=matrix(0,nrow=dim(A3)[1],ncol=R)
    for(r in 1:R){
Grad[,r]=ttl(as.tensor(tem),list(as.matrix(t(A1[,r])),as.matrix(t(A2[,r]))),ms=c(1,2))@data
    }}else if(mode==2){
    Grad=matrix(0,nrow=dim(A2)[1],ncol=R)
    for(r in 1:R){
Grad[,r]=ttl(as.tensor(tem),list(as.matrix(t(A1[,r])),as.matrix(t(A3[,r]))),ms=c(1,3))@data
    }}else if(mode==1){
    Grad=matrix(0,nrow=dim(A1)[1],ncol=R)
    for(r in 1:R){
Grad[,r]=ttl(as.tensor(tem),list(as.matrix(t(A2[,r])),as.matrix(t(A3[,r]))),ms=c(2,3))@data
            }}
     return(Grad)
}

cost=function(A1,A2,A3,Ybar,W,type=c("logistic","hinge")){
    return(mean(W*loss(tensorize(A1,A2,A3)*Ybar,type),na.rm=TRUE))
}
loss=function(y,type=c("logistic","hinge")){
    if(type=="hinge") return(ifelse(1-y>0,1-y,0))
    if(type=="logistic") return(log(1+exp(-y)))
}

### ADMM for classification
ADMM=function(Ybar,W,r,rho=0.1,lambda=10^(-3)){
  result=list();
  
  Lambda=array(0,dim=dim(Ybar))
  PQ=0; iter=0; obj =residual=error=max.thresh;
  
  rho_list=NULL;
  while(((iter < 10)|(error > 10^-3))){

    PQ_prev=PQ
    
    ### update B
    if((rho+lambda)!=0){
    res=SVM_offset(Ybar,W,OffsetC=(2*rho*PQ-Lambda)/(2*(lambda+rho)),cost=1/(2*(rho+lambda)))
    }else if((rho+lambda)==0){
         res=SVM_offset(Ybar,W,OffsetC=array(0,dim=dim(W)),cost=max.thresh)
    }
    obj=c(obj,res$hinge) ## minimize objective
    B=res$coef
    
    ## Update PQ
    if(rho==0){PQ=B}
    else{
    if(length(dim(Ybar))>2){
   sink("NULL")
    PQ=cp(as.tensor(B+1/(2*rho)*Lambda),r)
    sink()
        PQ=PQ$est@data
    }else if(length(dim(Ybar))==2){
        PQ=svd(B+1/(2*rho)*Lambda)
        if(r==1){
            PQ=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
        }else{
            PQ=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
        }
    }
    }
    
    
    residual=c(residual,sqrt(sum((B-PQ)^2)))
    
    ## update Lambda
    Lambda=Lambda+2*rho*(B-PQ)
    
    ## geometric step size
    if(iter>=10){
        rho=rho*1.1;
        lambda=lambda*1.1;
    }
    

    rho_list=c(rho_list,rho)
    iter=iter+1;

    error=abs(-residual[iter+1]+residual[iter])
    if(iter>=50) break
  }
  
  result$obj=obj[-1];
  result$iter=iter;
  result$error=error;
  result$fitted=PQ; ## exact low-rank
  result$B=B; ## approximate low-rank from SVM
  result$residual=residual[-1];result$rho=rho_list;
  return(result)
}

cost_svm=function(x,Ybar,W,OffsetC,cost){
   x=array(x,dim=dim(Ybar))
   return(mean(W*loss((x+OffsetC)*Ybar,"hinge"),na.rm=TRUE)+1/(2*cost)*mean(x^2,na.rm=TRUE))
}
#grad_svm=function(x,Ybar,W,OffsetC,cost){
#   x=array(x,dim=dim(Ybar))
#   margin=Ybar*(x+OffsetC)
#   tem=-W*Ybar*(margin<1)
#   return(tem+1/cost*x)
#}

SVM_offset=function(Ybar,W,OffsetC,cost=1,option=2){
    ### option 1
    if(option==1){## memory-saving option. did not check carefully
    coef=array(0,dim=dim(Ybar))
    opt=optim(coef,function(x)cost_svm(x,Ybar,W,OffsetC,cost))
#opt=optim(coef,function(x)cost_svm(x,Ybar,W,OffsetC,cost),function(x)grad_svm(x,Ybar,W,OffsetC,cost),method="BFGS")
coef=array(opt$par,dim=dim(Ybar))+OffsetC

return(list("res"=opt$value,"coef"=coef,"hinge"=objective(coef[nonmissing],Ybar[nonmissing],W[nonmissing])))
    }else{
### option 2
  n=length(Ybar)
  missing=which(is.na(Ybar)==TRUE)
  nonmissing=setdiff(1:n,missing)
  
  m=length(Ybar[nonmissing])
  dvec = 1-c(Ybar[nonmissing]*OffsetC[nonmissing])
  Dmat = diag(1,m)
  Amat = cbind(c(Ybar[nonmissing]),diag(1,m),-diag(1,m))
  bvec = c(rep(0,1+m),-c(cost*W[nonmissing]))
  res = solve.QP(Dmat,dvec,Amat,bvec,meq =1)
  
  ## calculate coefficient
  coef=OffsetC
  coef[nonmissing]=coef[nonmissing]+res$solution*Ybar[nonmissing]

  return(list("res"=res$value,"coef"=coef,"hinge"=objective(coef[nonmissing],Ybar[nonmissing],W[nonmissing])))
    }
}

hinge = function(y) ifelse(1-y>0,1-y,0)

objective=function(yfit,Ybar,W){
    return(sum(hinge(Ybar*yfit)*W))
}

likelihood = function(data,theta){
    index=which(is.na(data)==F & is.na(theta)==F)
   return(sqrt(sum((data[index]-theta[index])^2)))
}

##################### construct CP tensor using factor matrices X, Y, Z ###################################
tensorize=function(X,Y,Z){
    r=dim(X)[2]
    tensor=0
    if(is.matrix(X)==0){
        tensor=X%o%Y%o%Z
        return(tensor)
    }
    
    for(i in 1:r){
        tensor=tensor+X[,i]%o%Y[,i]%o%Z[,i]
    }
    return(tensor)
}

fit_continuous=function(data,r){
    index=which(is.na(data)==TRUE)
    data[index]=mean(data,na.rm=TRUE)
   original_data=data
   
    if(length(dim(data))>=3){
sink("NULL")
decomp=tryCatch(cp(as.tensor(original_data),r),error=function(c)"degeneracy")
suppressWarnings(sink())
if(inherits(decomp,"character")==TRUE){
    U=list();
    U[[1]]=matrix(0,ncol=r,nrow=dim(data)[1])
    U[[2]]=matrix(0,ncol=r,nrow=dim(data)[2])
    U[[3]]=matrix(0,ncol=r,nrow=dim(data)[3])
return(list("est"=array(NA,dim=dim(data)),"U"=U,"lambda"=rep(0,r),"info"="degeneracy"))
}
    res0=1
    res=0
    thresh=10^(-3)
    error=NULL
    
    while((res0-res)>thresh){
    res0=likelihood(original_data,decomp$est@data)
    sink("NULL")
    decomp=cp(as.tensor(data),r)
    sink()
    res=likelihood(original_data,decomp$est@data)
    data[index]=decomp$est@data[index]
    error=c(error,res)
    }

    return(list("est"=decomp$est@data,"U"=decomp$U,"lambda"=decomp$lambda,"info"=decomp))
    }else if(length(dim(data))==2){
        PQ=svd(data)
        if(r==1){
            decomp=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
        }else{
            decomp=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
        }
        res0=1
        res=0
        thresh=10^(-3)
        error=NULL
        
        while((res0-res)>thresh){
            res0=likelihood(original_data,decomp)
            if(r==1){
                decomp=cbind(PQ$u[,1:r])%*%diag(as.matrix(PQ$d[1:r]))%*%t(cbind(PQ$v[,1:r]))
            }else{
                decomp=PQ$u[,1:r]%*%diag(PQ$d[1:r])%*%t(PQ$v[,1:r])
            }
            res=likelihood(original_data,decomp)
            data[index]=decomp[index]
            error=c(error,res)
        }
        return(decomp)
}
}

#################### simulation model ####################
graphon_to_tensor=function(a,b,c,type){
    d1=length(a);d2=length(b);d3=length(c)
    M=array(0,dim=c(d1,d2,d3))
    if(type==10){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=log(0.5+max(a[i],b[j],c[k]))
                }
            }
        }
    }
    if(type==9){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=2-exp(min(a[i],b[j],c[k])^(1/3)) ##
                }
            }
        }
    }
    if(type==6){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=abs(a[i]-b[j]) ## full rank
                }
            }
        }
    }
    if(type==7){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=1/(1+exp(max(a[i],b[j],c[k])+min(a[i],b[j],c[k])))
                }
            }
        }
    }
    if(type==8){
        for(i in 1:d1){
            for(j in 1:d2){
                for(k in 1:d3){
                    M[i,j,k]=exp(-max(a[i],b[j],c[k])^(3/4))
                }
            }
        }
    }
    if(type==5){ ## stochastic block model
        r1=(1:length(a))%%3+1
        r2=(1:length(b))%%3+1
        r3=(1:length(c))%%3+1
        set.seed(1)
        value=array(rnorm(3^3,0,1),dim=rep(3,3))
        value=(value-min(value))/(max(value)-min(value))
        M=value[r1,r2,r3]
    }
    if(type==1){
        T=tensorize(a,b,c)
        M=1/(1+exp(-100*T)) ## single index
    }
    if(type==2){

    M=tensorize(a,b,c)
    #M=M*M*M
        ## low rank model
    }
    return(M)
}
appx_rank=function(tensor,thresh=95,step=5){
    size=dim(tensor)
    size=sort(size)
    min=which(sqrt(cumsum(svd(unfold(tensor,1:2,3)@data)$d^2))>thresh*0.01)[1]
    res=test=NULL
    r=min
    while(r<=(size[1]*size[2])){
        r=r+step;
        sink("NULL");
        rank=c(r,cp(tensor,r)$norm_percent)
        sink();
        res=rbind(res,rank)
        if(rank[2]>thresh) break
    }
    return(res)
}
