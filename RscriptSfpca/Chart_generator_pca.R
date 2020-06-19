rm(list = ls())

##### To be set
choose_Scenario = 3  ### 1 , 2  or 3, TO BE SET

### You can run the following
### Generate charts for analysis of both datasets
source("UTILITIES.R")
graphics.off()

load("PcaResults.rdata")
attach(PcaData)
data_to_attach = forward 
attach(data_to_attach) 

### Legend 
ncolors = 1 + ndefecttype ###per la legenda
defect_names = c("OOC-WS","OOC-MS","OOC-IRR")
img_dir_names = c("img_Scenario_I","img_Scenario_II","img_Scenario_III")
legend = c("in-control",defect_names[choose_Scenario]) ## Change names when needed
img_dir = img_dir_names[choose_Scenario]
colors = c("gray","blue")
colvec = rep(colors[1],ncontrol)
lwds = rep(1,ncontrol)
for ( j in 1:ndefecttype){
   colvec = c(colvec,rep(colors[j+1],defectnumerosity[j])) 
   lwds = c(lwds, rep(1,defectnumerosity[j]))
}


#### Plot of the smoothed density functions
setwd(img_dir)
pdf(file = "PDFs_fP.pdf")
par(cex = 1.5)
matplot(t, smooth_dens, type='l', ylim=c(0,0.6), col=colvec, lty=1,ylab = "Bernstein Density Functions",xlab = "log(distance)",lwd = lwds)
legend("topleft",legend = legend,fill = colors, cex=0.7)
dev.off()
### Plot of projections and comparison
pdf(file = "PDFs_and_projections_fP.pdf",width=14, heigh=7)
par(mfrow=c(1,2), cex = 2)
matplot(t, densp, type='l', ylim=c(0,0.6), col=colvec, lty=1,ylab = paste0("Projected Density Functions k =",k),xlab = "log(distance)",lwd = lwds)
legend(x = a, y = 0.6,legend = legend,fill = colors, cex=0.7)
matplot(t, smooth_dens, type='l', ylim=c(0,0.6), col=colvec, lty=1,ylab = "Bernstein Density Functions",xlab = "log(distance)",lwd = lwds)
legend(x = a, y = 0.6,legend = legend,fill = colors, cex=0.7)
dev.off()

###T2 chart
pdf(file = "T2_fP.pdf")
par(cex = 1.5)
plot(1:ntrain,T2, pch=19,cex.lab = 1, cex.axis = 1, type='b',col = "black", xlim = c(1,ntrain), xlab = "element index", main = "T2 statistic", ylim = c(0,20),lwd = 2)
points(1:ntrain,T2,col = colvec, cex = 1.4,pch = 16)
abline(h = T2lim, col = "black")
legend("topleft",legend = legend,fill = colors, cex=0.7)
dev.off()
###Q chart 
pdf(file = "Q_fP.pdf")
par(cex = 1.5)
plot(1:ntrain,Q, pch=19,cex.lab = 1, cex.axis = 1, type='b',col = "black", xlim = c(1,ntrain), xlab = "element index", main = "Q statistic", ylim = c(0,100),lwd = 2)
points(1:ntrain,Q,col = colvec, cex = 1.4,pch = 16)
abline(h = Qlim, col = "black")
legend("topleft",legend = legend,fill = colors, cex=0.7)
dev.off()

#### Eigenfunctions plot
pdf(file = "PC1_fP.pdf")
par(cex = 1.5)
t_step = t[2]-t[1]
par(cex = 1.2)
plot(t,clr2density(lmean,t,t_step),type = "l",lwd = 4,ylim = c(0,0.7),ylab = "Perturbations",xlab = "log(distance)",main = paste("PC",1))
lines(t,clr2density(lmean+eigenf[,1]*3*sqrt(eigenv[1]),t,t_step),col = "red",lwd = 2)
lines(t,clr2density(lmean-eigenf[,1]*3*sqrt(eigenv[1]),t,t_step),col = "blue",lwd = 2)
dev.off()


detach(data_to_attach)



### REPEAT FOR fS

data_to_attach = backward
attach(data_to_attach)
pdf(file = "PDFs_fS.pdf")
par(cex = 1.5)
matplot(t, smooth_dens, type='l', ylim=c(0,0.6), col=colvec, lty=1,ylab = "Bernstein Density Functions",xlab = "log(distance)",lwd = lwds)
legend("topleft",legend = legend,fill = colors, cex=0.7)
dev.off()
### Plot of projections and comparison
pdf(file = "PDFs_and_projections_fS.pdf",width=14, heigh=7)
par(mfrow=c(1,2), cex = 2)
matplot(t, densp, type='l', ylim=c(0,0.6), col=colvec, lty=1,ylab = paste0("Projected Density Functions k =",k),xlab = "log(distance)",lwd = lwds)
legend(x = a, y = 0.6,legend = legend,fill = colors, cex=0.7)
matplot(t, smooth_dens, type='l', ylim=c(0,0.6), col=colvec, lty=1,ylab = "Bernstein Density Functions",xlab = "log(distance)",lwd = lwds)
legend(x = a, y = 0.6,legend = legend,fill = colors, cex=0.7)
dev.off()

###T2 chart
pdf(file = "T2_fS.pdf")
par(cex = 1.5)
plot(1:ntrain,T2, pch=19,cex.lab = 1, cex.axis = 1, type='b',col = "black", xlim = c(1,ntrain), xlab = "element index", main = "T2 statistic", ylim = c(0,20),lwd = 2)
points(1:ntrain,T2,col = colvec, cex = 1.4,pch = 16)
abline(h = T2lim, col = "black")
legend("topleft",legend = legend,fill = colors, cex=0.7)
dev.off()
###Q chart 
pdf(file = "Q_fS.pdf")
par(cex = 1.5)
plot(1:ntrain,Q, pch=19,cex.lab = 1, cex.axis = 1, type='b',col = "black", xlim = c(1,ntrain), xlab = "element index", main = "Q statistic", ylim = c(0,100),lwd = 2)
points(1:ntrain,Q,col = colvec, cex = 1.4,pch = 16)
abline(h = Qlim, col = "black")
legend("topleft",legend = legend,fill = colors, cex=0.7)
dev.off()

#### Eigenfunctions plot
pdf(file = "PC1_fS.pdf")
par(cex = 1.5)
t_step = t[2]-t[1]
par(cex = 1.2)
plot(t,clr2density(lmean,t,t_step),type = "l",lwd = 4,ylim = c(0,0.7),ylab = "Perturbations",xlab = "log(distance)",main = paste("PC",1))
lines(t,clr2density(lmean+eigenf[,1]*3*sqrt(eigenv[1]),t,t_step),col = "red",lwd = 2)
lines(t,clr2density(lmean-eigenf[,1]*3*sqrt(eigenv[1]),t,t_step),col = "blue",lwd = 2)
dev.off()

setwd("..")
detach(data_to_attach)
detach(PcaData)
