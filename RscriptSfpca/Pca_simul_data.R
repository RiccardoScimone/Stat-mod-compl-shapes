#### This script compute SFPCA for a dataset as in the paper. No changes needed
rm(list=ls())
require(foreach)
require(doParallel)
load("alldistances.rdata") ### Raw dataset

source("UTILITIES.r") ### Utilities for SFPCA

#### To reduce the numerosity of point clouds, average distances
mean_batch = function(v,batch){
  ret = rep(0,ceiling(length(v)/batch))
  for (i in 1:length(ret))
    ret[i]= mean(v[((i-1)*batch+1):min((i*batch),length(v))]);
  return(ret);
}

ncontrol = 95; ## number of control elements
ndefecttype = 1; ## number of simulated out of control elements
defectnumerosity = c(5);## number of defective elements per type
cutleft = 0.00002 #left cut for distance acquisition
cutright = 0.065#right cut for distance acquisition
fin_num = 500   ## number of distance samples after distance averaging (this is necessary to compute the bernstein polinomials)
thr = 0.98; #explained variance treshold

#### Set I type error probability
alfatot = 0.01
alfaside = 0.005
alfasingle = 1-sqrt(1-alfaside) ### alfa for a single control chart


ntrain = ncontrol + sum(defectnumerosity) ### total samples
a = log(cutleft)
b = log(cutright) ### a e b sono gli estremi dell'intervallo


data = list(); #### Prepare backward distance data for acquisition (dataset fS)

for (i in 1:ntrain){
  data[[i]] = alldistances[[i]]$backward
  data[[i]] = data[[i]][which(data[[i]]>cutleft & data[[i]]<cutright)]
}
data = lapply(data,sort) ### sort all distances
data = lapply(data,log) ### compute log of distances
batches = floor(as.numeric(lapply(data, length))/fin_num)
for (i in 1:ntrain)
  data[[i]] = mean_batch(data[[i]],batches[i]) ### distance averaging

#### Now we are ready for density estimation
t = seq(a,b,length=302)
t = t[-c(1,302)] ### support of all densities
ecdf = lapply(data,ecdf) ## COmpute empirical distribution function

### Prepare parallel environment
cl = makeCluster(6)
registerDoParallel(cl)

### Bernsein density estimators
epdf.sm1 = foreach(i = 1:ntrain)%dopar%{return(fd.01(ecdf[[i]],t,a,b))}
stopCluster(cl)


epdf.sm1a=simplify2array(epdf.sm1) ### from list to matrix for easier handling


### Now we perform SFPCA
t_step = t[2]-t[1] #### Step of the discretized interval
n=length(t) ### number of points
lepdf1 = lapply(epdf.sm1, clr, t, t_step) ###clr transform of densities
lepdf1a=simplify2array(lepdf1) ### ... now in matrix form

#### All pca is computed on clrs
m1=apply(lepdf1a, 1, mean) ### clr mean

#### we have to make some adjustement due to the trapezoidal formula
#### Used to approximate the L^2 scalar product
lmod=lepdf1a
lmod=lmod*sqrt(t_step)
lmod[1,]=lmod[1,]/sqrt(2)
lmod[n,]=lmod[n,]/sqrt(2)
S=(ntrain-1)/ntrain*cov(t(lmod)) ### Covariance matrix of clrs
pc=eigen(S) ### SVD of S 

### re-adjust principal components
pc$vec[1,]=pc$vec[1,]*sqrt(2)
pc$vec[n,]=pc$vec[n,]*sqrt(2)
pc$vec=pc$vec/sqrt(t_step)


Nmax.harm=10 #### Maximum number of eigenfunctions to consider

lepdf1a.c=NULL
for(i in 1:ntrain)
  lepdf1a.c=cbind(lepdf1a.c, lepdf1a[,i]-m1) #### centered clrs

# Compute the scores of the sample along the eigenfunctions
sc.a=matrix(NA, ncol=Nmax.harm, nrow=ntrain)
for(i in 1:ntrain)
{
  for(j in 1:Nmax.harm)
    sc.a[i,j]=trapzc(t_step, lepdf1a.c[,i]*pc$vec[,j])
}


#### Select the number of principal components
K=min(which((cumsum(pc$val)[1:Nmax.harm]/sum(pc$val)>thr))) 

### Compute the projections of the clrs on the selected principal subspace
lepdf1.p=NULL
for(i in 1:ntrain)
{
  tmp=m1
  for(j in 1:K)
    tmp=tmp+sc.a[i,j]*(pc$vec[,j])
  lepdf1.p=cbind(lepdf1.p, (tmp))
}

# map the projections to projections of densities
epdf1.p=NULL
for(i in 1:ntrain)
{
  epdf1.p=cbind(epdf1.p,clr2density(lepdf1.p[,i],t,t_step))
}

##### Now assemble the results
pca=list()
pca$values=pc$values ## eigenvalues
pca$harmonics=pc$vec ## eigenfunctions
pca$scores=sc.a
T2 = matrix(0, 1, ntrain) ### T2 statistics
Q =   matrix(0, 1, ntrain) ### Q statistics
SCO = rbind(sc.a) ### scores

### Compute the statistics
for (i in 1:ntrain)
{
  for (k in 1:K)
    T2[i] = T2[i]+(SCO[i,k]^2)/pc$values[k]
  for (k in (K+1):Nmax.harm)
    Q[i] = Q[i] + (SCO[i,k]^2)^2
}
T2 = as.vector(T2)

####Assemble control limits
T2lim = (ntrain-1)^2/ntrain*qbeta(1-alfasingle,K/2,(ntrain-K-1)/2)
### The limit for Q is more involved
th1 = sum(pca$values[(K+1):Nmax.harm]);
th2 = sum(pca$values[(K+1):Nmax.harm]^2)
th3 = sum(pca$values[(K+1):Nmax.harm]^3)
h0 = 1-2*th1*th3/(3*th2^2);
Qlim = th1*(1+qnorm(1-alfasingle)*sqrt(2*th2*h0^2)/th1+th2*h0*(h0-1)/th1^2)^(1/h0);

### put all in a list (it contains the result of SFPCA on the dataset fS)
backward = list(data,epdf.sm1a,lepdf1a,lepdf1.p,epdf1.p,m1,SCO,T2,Q,T2lim,Qlim,K,pca$harmonics,pca$values,t,n)
names(backward) = c("data","smooth_dens","ldens","ldensp","densp","lmean","scores","T2","Q","T2lim","Qlim","k","eigenf","eigenv","t","nt")


### Now we repeat exactly the same things for the dataset fp


data = list(); #### Prepare forward distance data for acquisition (dataset fP)

for (i in 1:ntrain){
  data[[i]] = alldistances[[i]]$forward
  data[[i]] = data[[i]][which(data[[i]]>cutleft & data[[i]]<cutright)]
}
data = lapply(data,sort) ### sort all distances
data = lapply(data,log) ### compute log of distances
batches = floor(as.numeric(lapply(data, length))/fin_num)
for (i in 1:ntrain)
  data[[i]] = mean_batch(data[[i]],batches[i]) ### distance averaging

#### Now we are ready for density estimation
t = seq(a,b,length=302)
t = t[-c(1,302)] ### support of all densities
ecdf = lapply(data,ecdf) ## COmpute empirical distribution function

### Prepare parallel environment
cl = makeCluster(6)
registerDoParallel(cl)

### Bernsein density estimators
epdf.sm1 = foreach(i = 1:ntrain)%dopar%{return(fd.01(ecdf[[i]],t,a,b))}
stopCluster(cl)


epdf.sm1a=simplify2array(epdf.sm1) ### from list to matrix for easier handling


### Now we perform SFPCA
t_step = t[2]-t[1] #### Step of the discretized interval
n=length(t) ### number of points
lepdf1 = lapply(epdf.sm1, clr, t, t_step) ###clr transform of densities
lepdf1a=simplify2array(lepdf1) ### ... now in matrix form

#### All pca is computed on clrs
m1=apply(lepdf1a, 1, mean) ### clr mean

#### we have to make some adjustement due to the trapezoidal formula
#### Used to approximate the L^2 scalar product
lmod=lepdf1a
lmod=lmod*sqrt(t_step)
lmod[1,]=lmod[1,]/sqrt(2)
lmod[n,]=lmod[n,]/sqrt(2)
S=(ntrain-1)/ntrain*cov(t(lmod)) ### Covariance matrix of clrs
pc=eigen(S) ### SVD of S 

### re-adjust principal components
pc$vec[1,]=pc$vec[1,]*sqrt(2)
pc$vec[n,]=pc$vec[n,]*sqrt(2)
pc$vec=pc$vec/sqrt(t_step)


Nmax.harm=10 #### Maximum number of eigenfunctions to consider

lepdf1a.c=NULL
for(i in 1:ntrain)
  lepdf1a.c=cbind(lepdf1a.c, lepdf1a[,i]-m1) #### centered clrs

# Compute the scores of the sample along the eigenfunctions
sc.a=matrix(NA, ncol=Nmax.harm, nrow=ntrain)
for(i in 1:ntrain)
{
  for(j in 1:Nmax.harm)
    sc.a[i,j]=trapzc(t_step, lepdf1a.c[,i]*pc$vec[,j])
}


#### Select the number of principal components
K=min(which((cumsum(pc$val)[1:Nmax.harm]/sum(pc$val)>thr))) 

### Compute the projections of the clrs on the selected principal subspace
lepdf1.p=NULL
for(i in 1:ntrain)
{
  tmp=m1
  for(j in 1:K)
    tmp=tmp+sc.a[i,j]*(pc$vec[,j])
  lepdf1.p=cbind(lepdf1.p, (tmp))
}

# map the projections to projections of densities
epdf1.p=NULL
for(i in 1:ntrain)
{
  epdf1.p=cbind(epdf1.p,clr2density(lepdf1.p[,i],t,t_step))
}

##### Now assemble the results
pca=list()
pca$values=pc$values ## eigenvalues
pca$harmonics=pc$vec ## eigenfunctions
pca$scores=sc.a
T2 = matrix(0, 1, ntrain) ### T2 statistics
Q =   matrix(0, 1, ntrain) ### Q statistics
SCO = rbind(sc.a) ### scores

### Compute the statistics
for (i in 1:ntrain)
{
  for (k in 1:K)
    T2[i] = T2[i]+(SCO[i,k]^2)/pc$values[k]
  for (k in (K+1):Nmax.harm)
    Q[i] = Q[i] + (SCO[i,k]^2)^2
}
T2 = as.vector(T2)

####Assemble control limits
T2lim = (ntrain-1)^2/ntrain*qbeta(1-alfasingle,K/2,(ntrain-K-1)/2)
### The limit for Q is more involved
th1 = sum(pca$values[(K+1):Nmax.harm]);
th2 = sum(pca$values[(K+1):Nmax.harm]^2)
th3 = sum(pca$values[(K+1):Nmax.harm]^3)
h0 = 1-2*th1*th3/(3*th2^2);
Qlim = th1*(1+qnorm(1-alfasingle)*sqrt(2*th2*h0^2)/th1+th2*h0*(h0-1)/th1^2)^(1/h0);

### put all in a list (it contains the result of SFPCA on the dataset fP)
forward = list(data,epdf.sm1a,lepdf1a,lepdf1.p,epdf1.p,m1,SCO,T2,Q,T2lim,Qlim,K,pca$harmonics,pca$values,t,n)
names(forward) = c("data","smooth_dens","ldens","ldensp","densp","lmean","scores","T2","Q","T2lim","Qlim","k","eigenf","eigenv","t","nt")





#### Put both results in a unique list, ready to be used for control
PcaData = list(ncontrol,ndefecttype,defectnumerosity,ntrain,a,b,thr,forward,backward)
names(PcaData) = c("ncontrol","ndefecttype","defectnumerosity","ntrain","a","b","thr","forward","backward")
save(PcaData,file = "PcaResults.rdata")

