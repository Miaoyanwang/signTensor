
library(imager)


file <- system.file('imgsample.jpg',package='imager')

im <- load.image("imgsample2.jpg")

plot(im,axes = F)
imarray = as.array(im)
dim(im)


### noise image
noise_img = imarray + array(rnorm(length(imarray),0,.1),dim = dim(imarray))
imnoise = as.cimg(noise_img)
plot(imnoise,axes = F)
denoise_img1 = denoise_img2 = noise_img


## missing image
missing=array(rbinom(length(imarray),1,1/2),dim=dim(imarray))
missing_img = imarray
missing_img[missing==1]=NA
immissing = as.cimg(missing_img)
plot(immissing,axes = F)
complete_img1 = complete_img2 = missing_img

## missing+noise image
par(mfrow = c(1,1))
mn_img = noise_img
mn_img[missing==1] = NA
immn = as.cimg(mn_img)
plot(immn,axes = F)
recover_img1 = recover_img2 = mn_img


########### noise image recovery####################################################
####################################################################################

noise_v = noise_img[,,1,]
str(noise_v)
Lmin = min(noise_v)
Lmax = max(noise_v)

################### method 1 ###################### 

noise_res = SignT(noise_v,10,Lmin=Lmin,Lmax=Lmax,H=10,option=3)
denoise_img1[,,1,] = noise_res$est
denoise_img1  = as.cimg(denoise_img1)
plot(denoise_img1,axes = F)


################### method 2 #######################
est2 = fit_continuous(noise_v,10)
denoise_img2[,,1,] = est2$est
denoise_img2  = as.cimg(denoise_img2)
plot(denoise_img2)


########### missing image recovery##################################################
####################################################################################


missing_v = missing_img[,,1,]
str(missing_v)
Lmin = min(missing_v,na.rm = T)
Lmax = max(missing_v,na.rm = T)

################## method 1 ########################

missing_res = SignT(missing_v,10,Lmin=Lmin,Lmax=Lmax,H=10,option=3)

complete_img1[,,1,] = missing_res$est
complete_img1  = as.cimg(complete_img1)
plot(complete_img1,axes = F)

plot(im,axes = F)
plot(immissing,axes = F)
################## method 2 ########################
est3 = fit_continuous(missing_v,10)
complete_img2[,,1,] = est3$est
complete_img2  = as.cimg(complete_img2)
plot(complete_img2,axes = F)





########### missing + noise image recovery##########################################
####################################################################################

mn_v = mn_img[,,1,]
str(mn_v)
Lmin = min(mn_v,na.rm = T)
Lmax = max(mn_v,na.rm = T)

################## method 1 ########################

mn_res = SignT(mn_v,10,Lmin=Lmin,Lmax=Lmax,H=10,option=3)
mn_res2 = SignT(mn_v,10,Lmin=Lmin,Lmax=Lmax,H=10,option=2)

recover_img1[,,1,] = mn_res$est
recover_img1  = as.cimg(recover_img1)

par(mfrow = c(1,1))
plot(recover_img1,axes = F)

plot(im,axes = F)
plot(immissing,axes = F)
plot(immn,axes = F)
################## method 2 ########################
est3 = fit_continuous(mn_v,10)
recover_img2[,,1,] = est3$est
recover_img2  = as.cimg(recover_img2)
plot(recover_img2,axes = F)


