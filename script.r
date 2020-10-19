library(poibin)
source("sharedfuns.r")
source("funnewww.r")
#############################################################
#############################################################
##### Step 1: Reading Data#####
retire.dat=retire.read.data()
report.delay=retire.dat$delay.prob
retire.dat.ld=retire.data.to.ld(retire.dat)
ER<-85
sigma.ret<-2/3
retire.obj=list(retire.dat.ld=retire.dat.ld,report.delay=report.delay,retire.pars=c(log(ER/gamma(1+sigma.ret)),sigma.ret))
################################################################################


################################################################################
##########Step 2: Estimation####
##### 2.1 Optimization####
minus.log.likehood.retire.dat.new(dat=retire.obj,pars=c(11,0.65))
fit1=lifetime.mle(dat=retire.obj,
minusloglik=minus.log.likehood.retire.dat.new, starts=c(7,-1))

save(fit1,file="fit.out.10000");
load("fit.out.10000");fit1


##### 2.2 Contour Ploting ####
zz2=ContourObj(6, 9, 50, .2, .6, 50, minus.log.likehood.retire.dat.new.contour.fun)
save(zz2,file="zz2.out");
load("zz2.out")
save(contour.fig,file="contour.fig.wwe")
fig.paper(filename="contour.fig.wwe", type="ps",width=7,height=7)
contour(zz2$x,zz2$y,zz2$z,levels=seq(min(zz2$z),min(zz2$z)+10,,10))#contour plotting
dev.off()

contour(zz2$x,zz2$y,zz2$z,levels=seq(min(zz2$z),min(zz2$z)+2.995732,,2))#contour plotting
points(thetahat.star.mat[,1],exp(thetahat.star.mat[,2]),pch=16,cex=.1,col=2)
contour(zz2$x,zz2$y,zz2$z,levels=seq(min(zz2$z),min(zz2$z)+2.995732,,2),lwd=2,add=T,col=3)

contour(zz2$x,zz2$y,zz2$z,levels=seq(min(zz2$z),min(zz2$z)+2.995732,,2))#contour plotting
points(theta.star[,1],exp(theta.star[,2]),pch=16,cex=.1,col=2)
contour(zz2$x,zz2$y,zz2$z,levels=seq(min(zz2$z),min(zz2$z)+2.995732,,2),lwd=2,add=T,col=3)

#remove zero rows
thetahat.star.mat=thetahat.star.mat[thetahat.star.mat[,1]!=0,]

#starts
system.time(tmp.pi<-retire.pi(tt=c(1,10,50,100,200,300,400),thetahat.mat=thetahat.star.mat,mle.obj=fit1,alpha=c(.025,.05,.5,.95,.975)))


############################################################################
############################################################################

#####Step 3: PREDICTION##################
#####3.1 reading data after UL#######
#retire.dat=retire.read.data()
#report.delay=retire.dat$delay.prob
#retire.dat.ld=retire.data.to.ld(retire.dat)
#retire.dat.ld1=retire.dat.ld
#retire.dat.ld1[33:46,1]=0  #set t_j^L to zero, assume all failures before t_j^L were reported.
#retire.dat.ld1[,4]=0
#retire.dat.ld1<-retire.dat.ld1[-c(1:32),]
#retire.obj1=list(retire.dat.ld=retire.dat.ld1,report.delay=report.delay,retire.pars=c(log(ER/gamma(1+sigma.ret)),sigma.ret))

################################################################################
######3.2 Point Prediction########
#predict.retire.dat.new(dat=retire.obj,pars=c(fit1$coef[1],exp(fit1$coef[2])),s=1)$num
#pp<-predict.retire.dat.new(dat=retire.obj,pars=c(fit1$coef[1],exp(fit1$coef[2])),s=1)$p


######3.3 Plotting##########################
cumulative.fail.num<-cumulative.obj(mle.obj=fit1,tt=1:400)#(n=400,dat=retire.obj,pars=c((fit1$coef[1]),exp(fit1$coef[2])))
save(cumulative.fail.num,file="cumulative.fail.num.10000")
#load("cumulative.fail.num.10000"); #cumulative.fail.num

fig.paper(filename="cumulative.fig1.wwe", type="ps",width=7,height=7)
cumulative.fig(n=400,obj=cumulative.fail.num)## the cumulative plotting
dev.off()

######3.4 PI Prediction####################
# pi.num(n=10000,s=6*20,mu=c(fit1$coef[1],fit1$coef[2]),sigma=fit1$vcov,
#dat=retire.obj,pars=c((fit1$coef[1]),exp(fit1$coef[2])),wts=retire.dat.ld1[,"Units"])
#test

zz<-fit(n=10000,B=1,mu=fit1$coef,sigma=fit1$vcov, dat=retire.obj1,pars=c((fit1$coef[1]),exp(fit1$coef[2])))
save(zz,file="pi.fig.data.out10000")

fig.paper(filename="fig1", type="ps",width=7,height=7)
pi.fig(n1<-20,n2<-400,obj1<-zz,obj2<-cumulative.fail.num)
dev.off()


#load(eta);#eta
#load(beta);#beta
#load("rho.pred1");#rho.pred1
#load("nstar");#nstar
#load("vb1");#vb1
#load("qq");#qq
#load("pp");#pp
#load("pi.fig.data.out10000");zz

##################################################################




