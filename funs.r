###Step1: Reading Data########
########################################################################
retire.read.data=function()
{
 failures=read.table("paper.failures.txt",sep="",header=T)
 unfailed=read.table("paper.unfail.txt",sep="",header=T)
 delay.prob=read.table("paper.report.delay.txt",sep="",header=T)
 res=list(failures=failures,unfailed=unfailed,delay.prob=delay.prob)
 return(invisible(res))
 print(res)
}
######################################################################
retire.data.to.ld=function(dat)
{
  failures=dat$failures 
  unfailed=dat$unfailed
  nn.f=dim(failures)[1]  
  nn.u=dim(unfailed)[1]
  Failed=c(rep(1,nn.f),rep(0,nn.u)) 
  TimeStart=c(failures[,"FailTime"]-0.5,rep(0,nn.u))  
  TimeEnd=c(failures[,"FailTime"]+0.5,unfailed[,"AgeFreeze"])
  AgeStart=c(rep(0,nn.f),rep(0,nn.u))
  AgeFreeze=c(failures[,"AgeAtFreeze"],unfailed[,"AgeFreeze"])
  Units=c(rep(1,nn.f),unfailed[,"NotReported"])
  
  dat.ld=cbind(TimeStart=TimeStart,TimeEnd=TimeEnd, Failed=Failed, AgeStart=AgeStart,AgeFreeze=AgeFreeze,Units=Units) #Take a sequence of vector, matrix or data frames arguments and combine by columns or rows, respectively. These are generic functions with methods for other R classes
  dat.ld=as.data.frame(dat.ld)
  return(invisible(dat.ld))
}
###############################################################################
###############################################################################


####STEP 2: Estimate Parameters#######

####(1)Likelihood Function#######
minus.log.likehood.retire.dat.new=function(dat,pars)
{
    dat.ld=dat$retire.dat.ld
    report.delay=dat$report.delay
    muR=dat$retire.pars[1]
    sigmaR=dat$retire.pars[2]

    mu=pars[1]
    sigma=exp(pars[2])

   #print(mu)
   #print(sigma)

    dl=dim(report.delay)[1]
    dd=report.delay[,1]
    prob.dd=report.delay[,2]

    tmp.f=function(x)
    {
     xx=log(x)
     zz=(xx-mu)/sigma
     zz1=(xx-muR)/sigmaR
     res=(1/(x*sigma))*dsev(zz)*(1-psev(zz1))
     #print(x)
 
     #print(res)
     return(res) 
    }


    n1=32
    n2=14

    mat1=matrix(0,nrow=n1,ncol=dl)
    #print(mat1)

    mat2=matrix(0,nrow=n2,ncol=dl)
    #print(mat2)

    for(i in 1:n1)
    {
      for(j in 1:dl)
      {
        t.int=int.bound(dat.ld[i,1],dat.ld[i,2],dat.ld[i,4],dat.ld[i,5],j-1)
         #print(c(i,t.int))
         #print(c(mu,sigma))
        t.res=integrate(tmp.f,lower=t.int[1], upper=t.int[2])[[1]]#,rel.tol = .Machine$double.eps^0.5)[[1]]
        mat1[i,j]=as.numeric(t.res)
         #print(mat1)
      }
    }

    mat1=mat1%*%prob.dd    #print(mat1)
    pp1=rowSums(mat1)
    #print(pp1)


    for(i in 1:n2)
    {
      for(j in 1:dl)
      {
        t.int=int.bound(dat.ld[n1+i,1],dat.ld[n1+i,2],dat.ld[n1+i,4],dat.ld[n1+i,5],j-1)
        t.res=integrate(tmp.f,lower=t.int[1], upper=t.int[2])[[1]]
        mat2[i,j]=as.numeric(t.res)
        #print(mat2)
        print(c(i,t.int))
      }
    }
       mat2=sweep(mat2,2,prob.dd,"*")
    #print(mat2)

    pp2=1-rowSums(mat2)

    tmp1=c(pp1,pp2)
    #print(tmp1)
    log.p=log(tmp1)
    log.p=log.p*dat.ld[,"Units"]
    mll=-sum(log.p)
    return(mll)
}
################################################################################
thetahat.star.random.wts=function(mle.obj=fit1,B=10)
{
  dat=mle.obj$dat
  dat.ld=dat$retire.dat.ld
  dat.wts=dat.ld[,"Units"]

  tmp.fun=function(x)
  {
    res=sum(rexp(x))
    return(res)
  }

  res=matrix(0,nrow=B,ncol=2)

  for(i in 1:B)
  {
    print(i)
    wts=apply(as.matrix(dat.wts),1,tmp.fun)
    dat.ld[,"Units"]=wts
    dat$retire.dat.ld=dat.ld

    fit.try=try(lifetime.mle(dat=dat, minusloglik=minus.log.likehood.retire.dat.new, starts=mle.obj$coef),silent=T)
    try.flag=attr(fit.try,"class")=="try-error"

    if(length(try.flag)==0)
    {
     res[i,]=fit.try$coef
    }
  }

  return(res)
}
####################################################################################
####(2) Int.Bound####
int.bound=function(intl,intu,tl,tu,dd)
{
  tld=intl+dd
  tud=intu+dd
  if((tud<=tl)|(tld>=tu))
  {
    res=c(intu,intu)
  }
  
  if((tud<=tu)&(tld>=tl))
  {
    res=c(intl,intu)
  }
  
  if(tud<=tu&tud>tl&tld<tl)
  {
    res=c(tl-dd,intu)
  }

  if(tud>tu&tld<tu&tld>=tl)
  {
    res=c(intl,tu-dd)
  }

  if(tud>tu&tld<tl)
  {
    res=c(tl-dd,tu-dd)
  }
  
  return(res)
}
###########################################################################
####(3)Optimization####
lifetime.mle=function(dat, minusloglik, starts, method = "BFGS")
{
  call=match.call()
  f = function(p) {
    minusloglik(dat,p) 
    }
  oout = optim(starts, f, method = method, hessian = TRUE)#,control=list(trace=T))
  coef = oout$par
  vcov =solve(oout$hessian)
  min = oout$value
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat))
}

################################################################
####Contour Ploting####
minus.log.likehood.retire.dat.new.contour.fun=function(mu,sigma)
{
   res=minus.log.likehood.retire.dat.new(dat=retire.obj,pars=c(mu,log(sigma)))
   return(res)
}
#################################################################################
#################################################################################


#############Step3: Prediction#############################
###########################################################################
#####(1)Point Prediction####
predict.retire.dat.new=function(dat,pars,s)
{
    dat.ld=dat$retire.dat.ld
    report.delay=dat$report.delay
    muR=dat$retire.pars[1]
    sigmaR=dat$retire.pars[2]

    mu=pars[1]
    sigma=pars[2]
	
    dl=dim(report.delay)[1]
    dd=report.delay[,1]
    prob.dd=report.delay[,2]

    tmp.f=function(x)
    {
     xx=log(x)
     zz=(xx-mu)/sigma
     zz1=(xx-muR)/sigmaR
     res=(1/(x*sigma))*dsev(zz)*(1-psev(zz1))
     #print(x)
 
     #print(res)
     return(res) 
    }

    n1=14
    mat1=matrix(0,nrow=n1,ncol=dl)
    mat2=matrix(0,nrow=n1,ncol=dl)
       
  
    for(i in 1:n1)
    {
      for(j in 1:dl)
      {
        t.int=int.bound(dat.ld[i,1],dat.ld[i,2],dat.ld[i,4],dat.ld[i,5],j-1)
        #print(c(dat.ld[i,1],dat.ld[i,2],dat.ld[i,4],dat.ld[i,5],j-1,t.int))
        #t.res=integrate(tmp.f,lower=t.int[2]+1-0.5-(j-1), upper=t.int[2]+s+0.5-(j-1))[[1]]
        t.res=integrate(tmp.f,lower=t.int[2]+0.5, upper=t.int[2]+s+0.5)[[1]]
        #print(c(t.int[2]+1-0.5-(j-1), t.int[2]+s+0.5-(j-1),t.int[2]+0.5, t.int[2]+s+0.5))
        #print(c(dat.ld[i,4],dat.ld[i,5],j-1,t.int[2]+0.5, t.int[2]+s+0.5))
        mat1[i,j]=as.numeric(t.res)
        #print(mat1)
        #print(c(i,t.int))
        #print(j)
      }
    }

    mat1=sweep(mat1,2,prob.dd,"*")
   

    for(i in 1:n1)
    {
      for(j in 1:dl)
      {
        t.int=int.bound(dat.ld[i,1],dat.ld[i,2],dat.ld[i,4],dat.ld[i,5],j-1)
        #print(c(dat.ld[i,1],dat.ld[i,2],dat.ld[i,4],dat.ld[i,5],j-1,t.int))
        t.res=integrate(tmp.f,lower=t.int[1], upper=t.int[2])[[1]]
        mat2[i,j]=as.numeric(t.res)
        #print(mat2)
        #print(c(i,t.int))
      }
    }
    mat2=sweep(mat2,2,prob.dd,"*")
 
    pp1=rowSums(mat1)
    pp2=1-rowSums(mat2)
    #print(pp1)
    #print(pp2)
    p=pp1/pp2
    p1=p*(dat.ld[,"Units"])
   
    #print(p)
    #print(p1)
    mll=sum(p1)
    outcome<-list(prob=p,num=mll)
    return(outcome)
}

#######################ploting the cumulative curlve######################################
cumulative.obj<-function(mle.obj=fit1,tt)
{
  #x<-rep(1:n)
  dat=mle.obj$dat
  dat$retire.dat.ld<-dat$retire.dat.ld[-c(1:32),]
  
  nn=length(tt)
  x=double(nn)
  for(i in 1:nn)
  {
    x[i]<-predict.retire.dat.new(dat=dat,pars=c(mle.obj$coef[1],exp(mle.obj$coef[2])),s=tt[i])$num
  }
  res=cbind(tt,x)
  return(res)
}
################################################################################
cumulative.fig<-function(n,obj)
{
  plot(x<-seq(0,(n-1),by=1),y<-obj,type="l",xlab="weeks after DFT", ylab="Cumulative Number of Failure Time")
  lines(x<-seq(0,(n-1),by=1),y<-obj,col=3,lwd=1,lty=1)
}
#####################PI Simulation##############################
################################################################################
#####################PI Simulation##############################
#######################ploting the cumulative curlve######################################
cumulative.obj<-function(n,dat,pars)
{
  x<-rep(1:n)
  for(i in 1:n)
  {
   x[i]<-predict.retire.dat.new(dat=dat,pars=pars,s=i)$num
  }
  return(x)
}
################################################################################
cumulative.fig<-function(n,obj)
{
  plot(x<-seq(1,n,by=1),y<-obj,type="l",xlab="weeks after DFT", ylab="Cumulative Number of Failure Time")
  lines(x<-seq(1,n,by=1),y<-obj,col=3,lwd=1,lty=1)
}
################################################################################
#####################PI Simulation##############################
################################################################################
pi.num<-function(n,s,mu,sigma,dat,pars,wts,alpha=.05)
{
 #n is the number of repeats, same as B in the notes
 theta.star<-multi.norm.sim(n,mu,sigma,scale=1)
 save(theta.star,file="theta.star")
 mu.star<-theta.star[,1]
 sigma.star<-exp(theta.star[,2])

 #eta<-exp(mu)
 #beta<-1/exp(logsigma)
 #save(eta,file="eta")
 #save(beta,file="beta")

 ###########

 rho.pred<-function(n,s)
 {
   mat<-matrix(0,nrow=n,ncol=14)
   for (i in 1:n)
   {
     #print(i)
     x<-predict.retire.dat.new(dat=dat,pars=c(mu.star[i],sigma.star[i]),s=s)$p
     #print(x)
     mat[i,]<-x#as.numeric(x)
   }
   return(mat)
 }
 
 rho.pred1<-rho.pred(n,s)
 save(rho.pred1,file="rho.pred1")

 ############


 pp<-predict.retire.dat.new(dat=dat,pars=pars,s=s)$p
 #nstar=rpoibin(m=1000, pp=pp,wts=wts)
 nstar=rpoibin(m=n, pp=pp,wts=wts)

 save(pp,file="pp")

 #print(nstar)
 save(nstar,file="nstar")

 ###############
 vb<-function(n)
 {
  x<-double(n)
  for(i in 1:n)
  {
   x[i]<-ppoibin(kk=nstar[i],pp=rho.pred1[i,],method = "DFT-CF",wts=wts)
  }
  #pp<-print(pp)
  #save(pp,file="pp")
  return(x)
 }

 vb1<-vb(n)
 #print(vb1)
 save(vb1,file="vb1")

 ###################
 qq<-quantile(vb1,prob=c(alpha/2,1-alpha/2))

 print(qq)
 #save(qq,file="qq")

 ####################
 pp<-qpoibin(c(qq[1],qq[2]),pp=pp,wts=wts)#retire.dat.ld1[,"Units"])
 #print(pp)
 #save(pp,file="pp")

 return(pp)
}
################################################################################
fit<-function(n,B,mu,sigma,dat,pars,alpha=.05)
{
 #n is the number of repeats in simulation
 #B is the number of points
 mat1=matrix(0,nrow=B,ncol=2)
 #x=matrix(0,nrow=B,ncol=2)
 #for (i in 1:B)
 wts=dat$retire.dat.ld1[,"Units"]
 for (i in B)
 {
  x<-pi.num(n=n,mu=mu,sigma=sigma,dat=dat,par=pars,wts=wts,s=i*20,alpha=alpha)
  mat1[i,]=as.numeric(x)
 }
 return(mat1)
}
################################################################################
retire.pi<-function(tt,thetahat.mat,mle.obj,alpha=c(.025,.05,.5,.95,.975))
{
 nn=length(tt)
 mat1=matrix(0,nrow=nn,ncol=length(alpha))
 #x=matrix(0,nrow=B,ncol=2)
 #for (i in 1:B)
 #wts=dat$retire.dat.ld1[,"Units"]
 for (i in 1:nn)
 {
  x<-retire.pi.tt(tt=tt[i],thetahat.mat=thetahat.mat,mle.obj=mle.obj,alpha=alpha)
  mat1[i,]=as.numeric(x)
 }
 return(mat1)
}
################################################################################
retire.pi.tt<-function(tt,thetahat.mat,mle.obj,alpha)
{
 #n,s,mu,sigma,dat,pars,wts,alpha=.05)
 mu.star<-thetahat.mat[,1]
 sigma.star<-exp(thetahat.mat[,2])
 dat=mle.obj$dat
 dat$retire.dat.ld=dat$retire.dat.ld[-c(1:32),]
 xpars=mle.obj$coef
 pars=c(xpars[1],exp(xpars[2]))
 wts=dat$retire.dat.ld[,"Units"]
 n=dim(thetahat.mat)[1]

 rho.pred1<-matrix(0,nrow=n,ncol=14)
 for (i in 1:n)
 {
  #print(i)
  x<-predict.retire.dat.new(dat=dat,pars=c(mu.star[i],sigma.star[i]),s=tt)$prob
  #print(x)
  rho.pred1[i,]<-x#as.numeric(x)
 }

 #save(rho.pred1,file="rho.pred1")
 ############


 pp<-predict.retire.dat.new(dat=dat,pars=pars,s=tt)$prob
 #nstar=rpoibin(m=1000, pp=pp,wts=wts)
 nstar=rpoibin(m=n, pp=pp,wts=wts)

 #save(pp,file="pp")

 #print(nstar)
 #save(nstar,file="nstar")

 ###############

 vb1<-double(n)
 for(i in 1:n)
 {
   vb1[i]<-ppoibin(kk=nstar[i],pp=rho.pred1[i,],method = "DFT-CF",wts=wts)
 }

 #save(vb1,file="vb1")

 ###################
 qq<-quantile(vb1,prob=alpha)

 print(qq)
 #save(qq,file="qq")

 ####################
 pp<-qpoibin(qq,pp=pp,wts=wts)#retire.dat.ld1[,"Units"])
 #print(pp)
 #save(pp,file="pp")

 return(pp)
}
################################################################################
pi.fig<-function(n1,n2,obj1,obj2)
{
 plot(x<-seq(1,20*(n1-1)+1,by=20),y<-obj1[,2],ylim=c(0,120),type="l",xlab="weeks after DFT", ylab="Cumulative Number of Failure Time")
 lines(x<-seq(1,20*(n1-1)+1,by=20),y<-obj1[,2],col=4,lty=4,lwd=2)
 #plot(x<-seq(1,20*(n1-1)+1,by=20),y<-obj1[,1],type="l",xlab="weeks after DFT", ylab="Cumulative Number of Failure Time")
 lines(x<-seq(1,20*(n1-1)+1,by=20),y<-obj1[,1],col=1,lty=4,lwd=2)
 #plot(x<-seq(1,n2,by=1),y<-obj2,type="l",xlab="weeks after DFT", ylab="Cumulative Number of Failure Time")
 lines(x<-seq(1,n2,by=1),y<-obj2,col=3,lwd=1,lty=1)
}
################################################################################


