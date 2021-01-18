source("signT.R")
library(rgl)
library("plot3D")

library("RColorBrewer")
marker = list(color = brewer.pal(11, "RdBu"))$color
d=9
par(mfrow = c(1, 2))
x <- y <- seq(0, 1, length=d)
z <- outer(x, y, function(x,y)log(1+0.5*pmax(x,y)))

pdf("image3.pdf",width=30,height=10)
set.seed(1)
par(mfrow=c(1,3))
data=z+rnorm(length(z),0,0.01)
missing=array(rbinom(length(z),1,0.05),dim=c(d,d,1))
data[missing==1]=NA
persp3D(z = data,col= marker,zlim=c(min(data,na.rm=T),max(data,na.rm=T)+0.02),border="gray",phi = 30,theta=10)

res=SignT(array(data,dim=c(d,d,1)),2,Lmin=min(z),Lmax=max(z),H=5,option=2)

image3D(z = 0, colvar = sign(res$fitted[,,1,2]),zlim=c(-3,20),col=c("gray","white"), border="black",alpha=0.7,theta=10)
image3D(z = 5, colvar = sign(res$fitted[,,1,4]),add=TRUE,col=c("gray","white"),border="black")
image3D(z = 10, colvar = sign(res$fitted[,,1,6]),add=TRUE,col=c("gray","white"),border="black")
image3D(z = 15, colvar = sign(res$fitted[,,1,8]),add=TRUE,col=c("gray","white"),border="black")
image(sign(res$fitted[,,1,8]))


image3D(z=0.5,colvar=z,col=marker,zlim=c(0,1),phi = 25,theta=30)
dev.off()

image(data,col=marker)
image(res$est[,,1],col=marker)
