# Demonstration code
# Paper: Inverse Models for Spatio-Temporal Processes using Spatially Distributed Sensor Data Streams
# Author: Xiao Liu
# Contact: XL027@uark.edu
# Webpage: https://sites.google.com/site/liuxiaosite1/home

# Description:
# This code allows users to: 
# 1). choose one of the monitoring network designs: irregular, non-uniform and shifted uniform
# 2). obtain sensor monitoring data
# 3). perform inverse modeling
# 4). visulize the results

# About the ADMM algorithm:
# 1). the code provides the ADMM algorithm using the regularizations described in the paper
# 2). the code also allows users to use existing "Fused Lasso" and "Elastic Net" from the BB package
#   [the BB package is for optimizing a high-dimensional nonlinear objective function]

# ------------------------------------
# ------------------------------------
# ------------------------------------
# Packages needed 
# ------------------------------------
# ------------------------------------
# ------------------------------------
rm(list=ls())
library(sp)
library(RColorBrewer)
library(spectral)
library(Matrix)
library(fpp2)
library(BB) 
# ------------------------------------
# ------------------------------------
# ------------------------------------
# Work directory and parameter settings
# ------------------------------------
# ------------------------------------
# ------------------------------------
wd = "/Users/Xiao/Desktop/code_inverse/demo_code_GitHub"
setwd(wd)
# some colors for visualization
colPalette.basis = adjustcolor(colorRampPalette(brewer.pal(n=9, 'Greys'))(100), .85)
colPalette = adjustcolor(colorRampPalette(brewer.pal(n=9, 'YlOrRd'))(100), .85)
colPalette2 = adjustcolor(colorRampPalette(brewer.pal(n=9, 'Blues'))(100), 0.99)
colPalette1 = rev(rainbow(34, start=5/6, end=3/6, alpha=0.8))
colPalette4 = colorRampPalette(c('blue', 'yellow','red', 'black'))(100)

# ------------------------------------
# ------------------------------------
# ------------------------------------
# Load ADMM sub-routines
# ------------------------------------
# ------------------------------------
# ------------------------------------
source("fused_L2_v6.R")

# ------------------------------------
# ------------------------------------
# ------------------------------------
# load data of the advection-diffusion process
# ------------------------------------
# ------------------------------------
# ------------------------------------
load("data_all_300_0.005_0.005_0_0.00025_1600.RData")
iii = 10 # only use data from t=1 to t=10 to perform the inverse modeling
N1 = N2 = sqrt(n.alpha)
N = N1*N2
k1 = rep( seq(-N1/2+1, N1/2, 1), each=N2)
k2 = rep( seq(-N2/2+1, N2/2, 1), N1)
data = data[,c(1:iii)]
L = ncol(data)
v = colMeans(mu)
D = sigma[[1]]
zeta = zeta.G
grd.x = unique(coordinates(grd.original)[,1])
grd.y = ( unique(coordinates(grd.original)[,2]) ) 

# plot the data loaded
t.display = 5 # show the underlying process at time t.display [t.display =1,2,...,iii]
data.spdf = SpatialPixelsDataFrame(grd.original, data.frame(data[,t.display ]))
plot(grd.original, axes=TRUE,col="grey",xlim=c(0,1),ylim=c(0,1))
image(data.spdf , add=TRUE, col=colPalette.basis)

# ------------------------------------
# ------------------------------------
# ------------------------------------
# create the sensor network
# the users may choose irregular, non-uniform or shifted uniform
# ------------------------------------
# ------------------------------------
# ------------------------------------

# If the users would like to generate irregular grid:
if (TRUE){
  M.s = 50 # specify the number of sensors in the network
  grd.candidate = cbind( rep( seq(0, 1-1/N1, 1/N1), each=N2), 
                         rep( seq(0, 1-1/N2, 1/N2), N1))
  set.seed(3)
  grd.sample = sample(c(1:nrow(grd.candidate)),M.s,replace = FALSE)
  grd.sample = sort(grd.sample)
  grd = grd.candidate[grd.sample,]
  grd.sp = SpatialPoints(grd)
  
  # plot the grid
  plot(grd.sp,axes=TRUE)
}

# If the users would like to generate the non-uniform grid, change FALSE to TRUE
if (FALSE){
  M.s = 10 # specify the number of sensors in a row and a column
  set.seed(3)
  tmp1 = sort( sample(seq(0, 1-1/40, 1/40), M.s, replace=FALSE) )
  tmp2 = sort( sample(seq(0, 1-1/40, 1/40), M.s, replace=FALSE) )
  grd.sample = cbind( rep( tmp1, each=M.s), 
                      rep( tmp2, M.s))
  grd = grd.sample
  grd.sp = SpatialPoints(grd)
  plot(grd.sp,axes=TRUE)
}

# If the users would like to generate the shifted-uniform grid, change FALSE to TRUE
if (FALSE){
  M.s =  5 # specify the dimension of a uniform grid
  grd1 = cbind( rep( seq(0, 1-1/M.s, 1/M.s), each=M.s), 
                rep( seq(0, 1-1/M.s, 1/M.s), M.s))
  grd1.sp = SpatialPoints(grd1)
  grd2 = cbind( rep( seq(0, 1-1/M.s, 1/M.s)+0.04, each=M.s), 
                rep( seq(0, 1-1/M.s, 1/M.s)+0.175, M.s)) 
  # the shift along x and y directions: 0.04 and 0.175
  grd2.sp = SpatialPoints(grd2)
  grd.sp = SpatialPoints(rbind(grd1, grd2))
  plot(grd.sp,axes=TRUE)
}

# ------------------------------------
# ------------------------------------
# ------------------------------------
# Obtain the sensor observations based ont the nework created above
# ------------------------------------
# ------------------------------------
# ------------------------------------
data.grd = array(0, dim=c(M.s,L))
for (l in 1:L){
  # data at grd1
  for (i in 1:length(grd.sp)){
    tmp1 = coordinates(grd.sp)[i,1] - coordinates(grd.original)[,1]
    tmp2 = coordinates(grd.sp)[i,2] - coordinates(grd.original)[,2]
    dist = sqrt( tmp1^2 + tmp2^2 )
    #if (min(dist)<= ((1/N1)^2+(1/N2)^2)/2){
    #  case = which(dist == min(dist))
    #  data.grd[i,l] = data[case,l]
    #}else{
    case = which(dist == min(dist))
    data.grd[i,l] = mean(data[case,l])
  }
  
  data.grd.spdf = SpatialPointsDataFrame(grd.sp, data.frame(data.grd[,l]))
}

# visualization: plot the underlying process (unobserved) and the observed data by selected sensors
par(mfcol=c(1,2))
par(mar=c(4,2,1,1))

# plot the underlying process (unobserved) 
sensor.select=c(10, 25, 35) # plot the sensor observations by sensors #10, #25, and #50
image(x = unique(coordinates(grd.original)[,1]),
      y = unique(coordinates(grd.original)[,2]),
      z = (t(matrix((data[,1]),nrow=n.y))),
      col=heat.colors(255), xlab="",ylab="")
arrows(x0= rep(seq(0.1,0.9,0.2),each=5), y0=rep(seq(0.1,0.9,0.2),5), 
       x1 =rep(seq(0.1,0.9,0.2),each=5)+0.05, y1 =rep(seq(0.1,0.9,0.2),5)+0.05,
       col="white",length=0.15,lwd=1)
plot(grd.sp,add=TRUE,pch=3,cex=0.5,col="darkgreen")
text(x=c(0.4,0.2,0.5),  y = c(0.2,0.4,0.5), labels = c(1:3),col="black")
text(x=coordinates(grd.sp)[ sensor.select,1]+0.02, 
     y = coordinates(grd.sp)[ sensor.select,2], 
     labels = c("A","B","C"),col="black")

# add the locations of the three selected sensors
plot(grd.sp[sensor.select],add=TRUE,pch=1,cex=2,col="blue")

# plot the sensor observation
set.seed(1)
plot(data.grd[sensor.select[1],]+rnorm(80,0,2),type="l",lwd=2,lty=1,
     ylim=c(0,200),xlab="time",ylab="observed values")
lines(data.grd[sensor.select[2],]+rnorm(80,0,2),lwd=2,lty=2)
lines(data.grd[sensor.select[3],]+rnorm(80,0,2),lwd=2,lty=3)
legend("topright",legend=paste("sensor: ",  c("A","B","C"), sep=""),
       lwd=rep(2,3),lty=c(1,2,3))


# ------------------------------------
# ------------------------------------
# ------------------------------------
# Inverse modeling based on the data observed
# ------------------------------------
# ------------------------------------
# ------------------------------------

# --------------------------
# Step 1: create the linear model
# --------------------------

# create Y
Y = matrix(data.grd, ncol=1)

# create X
q1 = rep( seq(-N1/2+1, N1/2, 1), each=N2)
q2 = rep( seq(-N2/2+1, N2/2, 1), N1)
q = cbind(q1, q2)
n.q = nrow(q)

# compute the G matrix
G.list = list()
gamma = array(0, dim=c(N, 1))
for (l in 1:L){
  for (i in 1:n.q){
    gamma[i,1] = complex(length.out = 0, real = - t(as.matrix(q[i,])) %*% D %*% as.matrix(q[i,])-zeta, 
                         imaginary = - t(as.matrix(v)) %*% as.matrix(q[i,]))
  }
  G.list[[l]] = exp( (l-1)*gamma*delta.t )
}

# compute the Fouier bases
FF = array(0/0, dim=c(M.s,N))
s = coordinates(grd.sp)
for (i in 1:M.s){
  for (j in 1:N){
    FF[i,j] = complex(length.out = 1, real = cos(2*pi*(s[i,1]*q[j,1]+s[i,2]*q[j,2])), 
                      imaginary = sin(2*pi*(s[i,1]*q[j,1]+s[i,2]*q[j,2])))
  }
}

# create the X matrix (design matrix)
X = array(0/0, dim=c(M.s*L, N))
for (i in 1:L){
  i.start = (i-1)*M.s + 1
  i.end = i*M.s 
  tmp = array(0, dim=c(N, N))
  for (j in 1:N){
    tmp[j,j] = G.list[[i]][j]
  }
  X[i.start:i.end, ] = FF %*% tmp
}

# --------------------------
# Step 2: create the regularizations
# --------------------------
w = 1
S.1 = w * diag(1,ncol(X)) 

# pairs in the x direction
pair.list = list()
j = 1
zero.vec = array(0,dim=c(1,n.q))
for (i in 1:(n.q-1)){
  tmp1 = q[i,1]
  tmp2 = q[i,2]
  q.target = which( (q[,1]==(tmp1+1))&(q[,2]==tmp2) )
  if (length(q.target)>0){
    tmp = zero.vec
    tmp[c(i, q.target)] = c(1,-1)
    pair.list[[j]] = tmp
    j = j + 1
  }
}
S.2 = do.call(rbind, pair.list)
q[which(pair.list[[200]]==-1),]

# pairs in the y direction
pair.list = list()
j = 1
zero.vec = array(0,dim=c(1,n.q))
for (i in 1:(n.q-1)){
  tmp1 = q[i,1]
  tmp2 = q[i,2]
  q.target = which( (q[,1]==(tmp1))&(q[,2]==(tmp2+1)) )
  if (length(q.target)>0){
    tmp = zero.vec
    tmp[c(i, q.target)] = c(1,-1)
    pair.list[[j]] = tmp
    j = j + 1
  }
}
S.3 = do.call(rbind, pair.list)

S = rbind(S.1, S.2, S.3)
SS = rbind(S.2, S.3)

# --------------------------
# Step 3: perform ADMM
# --------------------------
tmp1 = Re(X)
tmp2 = Im(X)
tmp.list = list()
x.ext = rbind( cbind(tmp1, -tmp2), cbind(tmp2, tmp1) )

y.ext = c(Re(Y),Im(Y))

tmp1 = Re(S)
tmp2 = Im(S)
tmp.list = list()
s.ext = rbind( cbind(tmp1, -tmp2), cbind(tmp2, tmp1) )

tmp1 = Re(SS)
tmp2 = Im(SS)
tmp.list = list()
ss.ext = rbind( cbind(tmp1, -tmp2), cbind(tmp2, tmp1) )

# Run the ADMM code written by the author (recommended)
# If the users do NOT want to run the ADMM code written by the authors, change TRUE to FALSE
if (TRUE){ # My ADMM
  cut.off = 4 # keep only the low frequency terms
  q.high.case = which( (abs(q[,1])<= cut.off)&(abs(q[,2])<= cut.off) )
  q.high.case = which(  sqrt( q[,1]^2 + q[,2]^2 ) <= cut.off )
  source("fused_L2_v6.R")
  tmp = fused.L2(y = y.ext, x = x.ext, s=ss.ext, r=1, 
                 lambda1=1, lambda2=1, q.high.case=q.high.case)
  
  eta.re = eta.im = array(0,dim=c(N,1))
  eta.re[q.high.case] = tmp[[1]][1:length(q.high.case)]
  eta.im[q.high.case] = tmp[[1]][(length(q.high.case)+1):(2*length(q.high.case))]
  
  eta = complex(length.out = N, real = eta.re, 
                imaginary = eta.im)
  eta.0 = eta
  
  Y.est.0 = X%*%eta.0
}

# If the users would like to use the existing Fused Lasso provided by R, change FALSE to TRUE
if (FALSE){ # Fused Lasso using ADMM package
  cut.off = 4
  q.high.case = which(  sqrt( q[,1]^2 + q[,2]^2 ) <= cut.off )
  
  tmp3 = admm.genlasso(
    A=x.ext[,c(q.high.case, q.high.case+ncol(x.ext)/2)],
    b=y.ext,
    D=s.ext[,c(q.high.case, q.high.case+ncol(s.ext)/2)],
    lambda = 1,
    #lambda2 = 50,
    rho = 1,
    abstol = 1e-04,
    reltol = 0.01,
    maxiter = 1000
  )
  
  eta.re = eta.im = array(0,dim=c(N,1))
  eta.re[q.high.case] = tmp3$x[1:(length(tmp3$x)/2)]
  eta.im[q.high.case] = tmp3$x[(length(tmp3$x)/2+1):(length(tmp3$x))]
  
  eta = complex(length.out = N, real = eta.re, 
                imaginary = eta.im)
  eta.1 = eta
  
  Y.est.1 = X%*%eta.1
}

# If the users would like to use the existing Elastic Net provided by R, change FALSE to TRUE
if (FALSE){ # Elastic Net using ADMM package
  cut.off = 4
  q.high.case = which(  sqrt( q[,1]^2 + q[,2]^2 ) <= cut.off )
  
  tmp3 = admm.enet(
    A=x.ext[,c(q.high.case, q.high.case+ncol(x.ext)/2)],
    b=y.ext,
    lambda1 = 1,
    lambda2 = 1,
    rho = 1,
    abstol = 1e-04,
    reltol = 0.01,
    maxiter = 1000
  )
  
  eta.re = eta.im = array(0,dim=c(N,1))
  eta.re[q.high.case] = tmp3$x[1:(length(tmp3$x)/2)]
  eta.im[q.high.case] = tmp3$x[(length(tmp3$x)/2+1):(length(tmp3$x))]
  
  eta = complex(length.out = N, real = eta.re, 
                imaginary = eta.im)
  eta.2 = eta
  
  Y.est.2 = X%*%eta.2
}

# --------------------------
# Step 4: Obtain the inverse model output 
# --------------------------
# create Fourier bases
FFF = array(0/0, dim=c(N,N))
ss = coordinates(grd.original)
for (i in 1:N){
  for (j in 1:N){
    FFF[i,j] = complex(length.out = 1, real = cos(2*pi*(ss[i,1]*q[j,1]+ss[i,2]*q[j,2])), 
                       imaginary = sin(2*pi*(ss[i,1]*q[j,1]+ss[i,2]*q[j,2])))
  }
}
image.recover.0 = FFF %*% eta.0
image.recover.m.0 = sqrt(Re(image.recover.0)^2 + Im(image.recover.0)^2)


# If the users used the Fused Lasso provide by R, change FALSE to TRUE
if (FALSE){
  image.recover.1 = FFF %*% eta.1
  image.recover.m.1 = sqrt(Re(image.recover.1)^2 + Im(image.recover.1)^2)
  image.recover.m.0 = image.recover.m.1
}
# If the users used the Elastic NEt provide by R, change FALSE to TRUE
if (FALSE){
  image.recover.2 = FFF %*% eta.2
  image.recover.m.2 = sqrt(Re(image.recover.2)^2 + Im(image.recover.2)^2)
  image.recover.m.0 = image.recover.m.2
}


# ------------------------------------
# ------------------------------------
# ------------------------------------
# Visualization of the output
# ------------------------------------
# ------------------------------------
# ------------------------------------

# plot the actual source center (unknown to the model)
source.center = cbind( c(0.4,0.2,0.5),c(0.2,0.4,0.5) )
source.center.sp = SpatialPoints(source.center)
plot(source.center, axes=TRUE, xlim=c(0,1), ylim=c(0,1), col="white",pch=20)
plot(grd.sp,add=TRUE,pch=3,cex=0.5,col="darkgreen")
# add contours of the model output
contour(x = unique(coordinates(grd.original)[,1]),
        y = unique(coordinates(grd.original)[,2]),
        z = (t(matrix((image.recover.m.0),nrow=n.y))),labels ="",
        add=TRUE,col="gray")
contour(x = unique(coordinates(grd.original)[,1]),
        y = unique(coordinates(grd.original)[,2]),
        z = (t(matrix((image.recover.m.0),nrow=n.y))),labels ="",
        levels = max(image.recover.m.0)*c(0.75, 0.85,0.95), add=TRUE,lwd=2,col="blue")
contour(x = unique(coordinates(grd.original)[,1]),
        y = unique(coordinates(grd.original)[,2]),
        z = (t(matrix((image.recover.m.0),nrow=n.y))),labels ="",
        levels = max(image.recover.m.0)*c(0.5,0.6), add=TRUE,lwd=2,col="orange",lty=2)
#plot(source.center.sp, add=TRUE, col="red",pch=20)
text(x=c(0.4,0.2,0.5),  y = c(0.2,0.4,0.5), labels = c(1:3),col="red")
title(paste("time: ", iii,sep=""))





