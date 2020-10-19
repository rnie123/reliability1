#created 08feb2010 for commonly used functions
################################################################################
#function for create a contour object # args of ff is not a vector
ContourObj=function(x1,x2,nx,y1,y2,ny,ff)
{
  f=function(tt)
  {
    x=tt[1]
    y=tt[2]
    return(ff(x,y))
  }
  xy.grid=list(x=seq(x1,x2,,nx),y=seq(y1,y2,,ny))
  t=expand.grid(xy.grid)
  zz=apply(as.matrix(t),1,f)
  z=matrix(zz,nx,ny,byrow=F)
  return(list(x=xy.grid$x,y=xy.grid$y,z=z))
}
################################################################################
#function for create a contour object, args of ff is a vector
ContourObj1=function(x1,x2,nx,y1,y2,ny,ff)
{
  xy.grid=list(x=seq(x1,x2,,nx),y=seq(y1,y2,,ny))
  tt=expand.grid(xy.grid)

  zz=apply(as.matrix(tt),1,ff)

  z=matrix(zz,nx,ny,byrow=F)
  return(list(x=xy.grid$x,y=xy.grid$y,z=z))
}
################################################################################
#output matrix to tex format
TexMat=function(x)
{
 M=as.matrix(x)
 m=nrow(M)
 n=ncol(M)
 cat("\\begin{pmatrix}",fill=T)
 for(i in 1:m)
 {
  tt=M[i,1]
  if (n>1)
  {
    for(j in 2:n)
    {
     tt=paste(tt,"&",M[i,j])
     }
   }
 cat(tt,"\\\\",fill=T)
 }
 cat("\\end{pmatrix}",fill=T)
}
################################################################################
TexMat1=function(x,nsmall)
{
 M=as.matrix(x)
 m=nrow(M)
 n=ncol(M)
 cat("\\begin{pmatrix}",fill=T)
 for(i in 1:m)
 {
  tt=format(M[i,1],nsmall=nsmall)
  if (n>1)
  {
    for(j in 2:n)
    {
     tt=paste(tt,"&",format(M[i,j],nsmall=nsmall))
     }
   }
 cat(tt,"\\\\",fill=T)
 }
 cat("\\end{pmatrix}",fill=T)
}

################################################################################
mat.unique=function(mat,cols)
{
 x=mat[,cols]
 temp <- apply(x, 1, function(x) paste(x, collapse = "\r"))
 dat=mat[!duplicated(as.vector(temp)),]
 return(dat)
}
#####cdf of sev#################################################################
psev=function(z)
{
  1-exp(-exp(z))
}
#####densigy of sev#############################################################
dsev=function(z)
{
  exp(z-exp(z))
}

################################################################################
################################################################################
dumpall=function()
{
  save(file="all.objs",list=ls(,envir=.GlobalEnv))
}
#####quantile of sev############################################################
qsev=function(p)
{
  p=ifelse(p>=.99999999999999,.99999999999999,p)
  p=ifelse(p<=1-.99999999999999,1-.99999999999999,p)
  log(-log(1-p))
}
################################################################################
#generate sev random variable
rsev=function(n)
{
  qsev(runif(n))
}
################################################################################
fig.paper=function(filename, type="ps",width=7,height=7)
{
  dir.create("./figures", showWarnings = F)
  if(type=="ps" | type=="eps")
  {
   postscript(file=paste("./figures/",filename,".",type,sep=""),width=width,
   height=height,horizontal=F)
  }
  if(type=="wmf")
  {
   win.metafile(filename = paste("./figures/",filename,".wmf",sep=""), 
   width = width, height = height, pointsize = 12, restoreConsole = TRUE)
  }

}
################################################################################
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
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat,minusloglik=minusloglik))
}
################################################################################
kaplan.meier.location=function(fit)
{
   xloc=fit$time[fit$n.event>0]
   yloc=1-fit$surv[fit$n.event>0]
   nn=length(yloc)
   yloc=(yloc+c(0,yloc[1:(nn-1)]))/2
   return(cbind(xloc,yloc))
}
################################################################################
miniusloglik.sev=function(dat, pars)
{
 #column 1 event time
 #column 2 for event indicator
 tt=dat[,1]
 delta=dat[,2]
 mu=pars[1]
 sigma=exp(pars[2])
 zz=(log(tt)-mu)/sigma
 ff=dsev(zz)/(sigma*tt)
 FF=psev(zz)
 
 ll=delta*log(ff)+(1-delta)*log(1-FF)
 res=(-1)*sum(ll)
 return(res)
}
################################################################################


################################################################################
##generate probability plot
################################################################################
SZ.probpaper=function(distribution="Weibull", xlab1 = "Time", ylab1="Fraction Failing", main1="Probability Paper",
x.range = c(1, 1000.), y.range = c(0.00001, 0.9) )
{
  probpaper(distribution=distribution, xlab = "", x.range = x.range, y.range =  y.range
        , original.par = F, my.title = NULL, cex = 1., cex.labs = 1, cex.tic.lab = 1, cex.title = cex.labs,
        sub.title = NULL, grids = 0., linear.axes = F, slope.axis = F, title.option = "paper", shape = NULL,
        draw.line = F, ylab = "",xlab1=xlab1,ylab1=ylab1,main1=main1)

}

probpaper=function(distribution, xlab = "", x.range = c(1, 1000.), y.range = c(0.000010999999999999999, 0.9
        ), original.par = F, my.title = NULL, cex = 1., cex.labs = 1, cex.tic.lab = 1, cex.title = cex.labs,
        sub.title = NULL, grids = 0., linear.axes = F, slope.axis = F, title.option = "paper", shape = NULL,
        draw.line = F, ylab = "",xlab1="Time",ylab1="Fraction Failing",main1="Probability Paper")
{
        #
        #
        #need to allow for sending down start values
        #
        old.par <- par()
        par(mar = c(5., 6., 4., 2.) + 0.10000000000000001, err = -1.)
        if(original.par)
                on.exit(par(old.par))
        if(is.null(my.title)) {
                my.title <- ""
        }
        xrna <- is.na(x.range)
        if(any(xrna))
                x.range[xrna] <- c(1., 10.)[xrna]
        yrna <- is.na(y.range)
        if(any(yrna))
                y.range[yrna] <- c(0.010999999999999999, 0.98999999999999999)[yrna]
        log.of.data <- probplot.setup(distribution, x.range, y.range, my.title = my.title, sub.title = "", cex =
                cex, ylab = ylab, xlab = "", grids = grids, linear.axes = linear.axes, title.option =
                title.option, draw.line = draw.line, slope.axis = slope.axis, shape = shape, cex.title =
                cex.title, cex.labs = cex.labs)
        title(main=main1,xlab = xlab, cex = cex.labs)
        mtext(ylab1, side = 2, line = -1,outer=T)
        mtext(xlab1, side = 1, line = -1.2,outer=T)

}


################################################################################
probplot.setup=function(distribution, x.range, y.range, xlab = "Time", cex = 1., cex.title = cex, shape = NULL, my.title = NULL,
        sub.title = "", grids = F, linear.axes = F, title.option = "blank", slope.axis = F, draw.line = F,
        cex.labs = 1., ylab = "", hw.xaxis = NULL, hw.yaxis = NULL, dump = T,
        cex.points = 1., title.line.adj = 0, ...)
{
        #
        #
        #set up for probability plotting
        #
        #cex.points is not really used here, but R gave warnings because
        #it was being passed in ... from the command six.npprobplot(ShockAbsorber.ld)
        #
        if(title.option != "blank" && (!is.logical(linear.axes) || (slope.axis))) {
                wqm.warning("Probably should not have a title if  shape parameter scale has been requested")
        }
        if(missing(x.range))
                wqm.warning("x.range should be specified")
        if(missing(y.range))
                wqm.warning("y.range should be specified")
        xdiff <- x.range[2.] - x.range[1.]
        if(xdiff < 0.)
                stop("Requested x-axis upper endpoint less than lower")
        x.range[1.] <- x.range[1.] + 9.9999999999999998e-013 * xdiff
        x.range[2.] <- x.range[2.] - 9.9999999999999998e-013 * xdiff
        ydiff <- y.range[2.] - y.range[1.]
        if(ydiff < 0.)
                stop("Requested y-axis upper endpoint less than lower")
        y.range[1.] <- y.range[1.] + 9.9999999999999998e-013 * ydiff
        y.range[2.] <- y.range[2.] - 9.9999999999999998e-013 * ydiff
        #
        #check and grab LockAxes, if requested
        #  also, store the current values for possible use later
        #
        GetAxesRange.out <- GetAxesRange("probplot.setup", x.axis = "xxx", x.range, xlab, y.axis = distribution,
                y.range, ylab)
        x.range <- GetAxesRange.out$x.range
        y.range <- GetAxesRange.out$y.range
        #
        #
        #this function sets par(mar=...). After this function is called,
        #there should be no more adjustments to par(mar) until there is
        #to be another high-level plot
        #
        #
        if(is.logical(linear.axes) && linear.axes) linear.axes <- "b"
        get.prob.scales.out <- get.prob.scales(distribution, shape = shape, prob.range = y.range)
        #       frame()
        log.of.data <- get.prob.scales.out$logger
        tick.label.loc <- as.numeric(get.prob.scales.out$tick.labels)
        tick.location <- as.numeric(get.prob.scales.out$tick.location)
        if(length(grep("Percent", ylab)) > 0)
                show.tick.labels <- get.prob.scales.out$percent.tick.labels
        else show.tick.labels <- get.prob.scales.out$tick.labels
        yp.range <- c(max(min(tick.label.loc), tick.label.loc[tick.label.loc < y.range[1.]]), min(max(
                tick.label.loc), tick.label.loc[tick.label.loc > y.range[2.]]))
        #probplot.setup.title.out <- probplot.setup.title(title.option = title.option, my.title = my.title,
        #        sub.title = sub.title, distribution = get.prob.scales.out$distribution, shape = shape)
        if((title.option == "paper" || title.option == "paper2") && grids != 0.) {
                #
                #
                #title option
                #   paper          xxx probability scale, cut top to allow scales
                #   paper2         no indication of which prob scale, cut top to allow scales
                #   full           full title (my.title, dist, and sub?)
                #   only.dist      only indicate the dist
                #   only.my        just my.title, omit dist
                #   only.sub       only subtitle
                #   blank          blank with no cut top
                #   blank2         blank, with cut.top
                #
                grids <- 2.
        }
        if(linear.axes == "q" || linear.axes == "b") {
                right.mar <- 6.
        }
        else {
                right.mar <- 2.
        }
        if(slope.axis) {
                left.mar <- 8.
        }
        else {
                left.mar <- 6.
        }
        old.par <- par(mar = c(3, left.mar-2, 3, right.mar) + 0.10000000000000001,
                err = -1.)
        y.range <- pp.quant(yp.range, distribution, shape)
        if(is.null(hw.xaxis)) {
                if(log.of.data)
                        data.axes.out <- logax(x.range[1.], x.range[2.], ...)
                else data.axes.out <- linax(x.range)
        }
        else {
                data.axes.out <- hw.xaxis
        }
        x.range <- pp.data(range(x.range, as.numeric(data.axes.out$ticlab)), log.of.data)
        #
        #set up the prob axes
        #
        if(is.null(hw.yaxis)) {
                tick.location <- tick.location[tick.location >= yp.range[1.] & tick.location <= yp.range[2.]]
                in.range <- tick.label.loc >= yp.range[1.] & tick.label.loc <= yp.range[2.]
                tick.labels <- get.prob.scales.out$tick.labels
        }
        else {
                tick.labels <- hw.yaxis$ticlab
                tick.label.loc <- as.numeric(hw.yaxis$ticlab)
                tick.location <- as.numeric(hw.yaxis$ticloc)
                tick.location <- tick.location[tick.location >= yp.range[1.] & tick.location <= yp.range[2.]]
                in.range <- tick.label.loc >= yp.range[1.] & tick.label.loc <= yp.range[2.]
        }
        #
        #make the plot outline
        #
        plot(x.range, y.range, xlab = "", ylab = "", xaxt = "n", yaxt = "n", type = "n", ...)
        title(xlab = xlab, cex = cex.labs)
        #
        #make the probability tick marks
        #
        axis(side = 2., at = pp.quant(tick.location, distribution, shape), labels = F, tck = -0.01, mgp = c(
                5., 1, 0.), cex = 1)
        #
        #make the probability tick labels
        #
        axis(side = 2., at = pp.quant(tick.label.loc[in.range], distribution, shape), labels = fix.exp.labels(
                show.tick.labels[in.range]), mgp = c(5., 0.7, 0.), cex =
                cex.labs, las = 1)
        #
        #set up linear axes corresponding to probability axes
        #
        if(linear.axes == "q" || linear.axes == "b") {
                lin.prob.axes.out <- linax(range(pp.quant(tick.location, distribution, shape)))
                data.tick.location <- as.numeric(lin.prob.axes.out$ticloc)
                data.tick.label.loc <- as.numeric(lin.prob.axes.out$ticlab)
                axis(side = 4., at = data.tick.location, labels = F, tck = -0.01, mgp = c(5., 2.1000000000000001,
                        0.), cex = 1.)
                axis(side = 4., at = data.tick.label.loc, labels = fix.exp.labels(lin.prob.axes.out$ticlab),
                        adj = 1., tck = -0.02, mgp = c(5., 2.3999999999999999, 0.), cex = cex.labs)
                mtext(side = 4., line = 3., text = paste("Standard", get.prob.scales.out$prob.scale, " quantile"),
                        cex = cex.labs)
        }
        data.tick.location <- as.numeric(data.axes.out$ticloc)
        data.tick.label.loc <- as.numeric(data.axes.out$ticlab)
        #
        #set up data tick marks
        #
        axis(side = 1., at = pp.data(data.tick.location, log.of.data), labels = F, tck = -0.01, mgp = c(5.,
                2.1, 0.), cex = 1.)
        xlabels <- vector.power10(data.axes.out$ticlab)
        #
        #set up data tick labels
        #
        #if(is.postsctiptok() && substring(xlabels[1.], 1., 1.) == "~") {
        #        mixed.mtext.vec(side = 1., at = pp.data(data.tick.label.loc, log.of.data), texts = xlabels, adj
        #                 = 0.5, cex = cex.labs, line = 0.75)
        #        axis(side = 1., at = pp.data(data.tick.label.loc, log.of.data), tck = -0.02, mgp = c(5., 1.,
        #                0.), labels = F, cex = cex.labs)
        #}
        #else {

                ##x axis
                #axis(side = 1., at = pp.data(data.tick.label.loc, log.of.data), labels = fix.exp.labels(
                #        data.axes.out$ticlab), adj = 0.5, tck = -0.02, mgp = c(5., 1., 0.), cex = cex.labs)
              aaass=fix.exp.labels(data.axes.out$ticlab)
              aaass=blank.fix(aaass)
              axis(side = 1., at = pp.data(data.tick.label.loc, log.of.data), 
              labels = aaass,mgp=c(3,0.7,0), cex = cex.labs)


        #}
        #
        #      do the grids
        #
        if(linear.axes == "t" || linear.axes == "b" && log.of.data) {
                lin.data.axes.out <- linax(range(log10(data.tick.location)))
                lin.data.tick.location <- as.numeric(lin.data.axes.out$ticloc) * 2.3025850000000001
                lin.data.tick.label.loc <- as.numeric(lin.data.axes.out$ticlab) * 2.3025850000000001
                axis(side = 3., at = lin.data.tick.location, labels = F, tck = -0.01, mgp = c(5.,
                        2.1000000000000001, 0.), cex = 1.)
                #       mtext(side = 3, line = 2.5, text = paste("Log", xlab), cex = 1.2)
                axis(side = 3., at = lin.data.tick.label.loc, labels = fix.exp.labels(lin.data.axes.out$ticlab),
                        adj = 0.5, tck = -0.02, mgp = c(5., 1., 0.), cex = cex.labs)
        }
        if(grids >= 1.) {
                #
                #vertical light grid lines
                #
                usr.out <- par("usr")
                yvec.low <- rep(usr.out[3.], length(data.tick.location))
                yvec.high <- rep(usr.out[4.], length(data.tick.location))
                #
                #horizontal light grid lines  tick.label.loc
                #
                #
                matlines(rbind(pp.data(data.tick.location, log.of.data), pp.data(data.tick.location, log.of.data)
                        ), rbind(yvec.low, yvec.high), col = 1., lty = 1., lwd = 0.20000000000000001)
                xvec.low <- rep(usr.out[1.], length(tick.location))
                xvec.high <- rep(usr.out[2.], length(tick.location))
                matlines(rbind(xvec.low, xvec.high), rbind(pp.quant(tick.location, distribution, shape), pp.quant(
                        tick.location, distribution, shape)), col = 1., lty = 1., lwd = 0.20000000000000001)
        }
        #
        #      write titles and axis lables
        #
        if(grids >= 2.) {
                #
                #vertical dark grid lines
                #
                yvec.low <- rep(usr.out[3.], length(data.tick.label.loc))
                yvec.high <- rep(usr.out[4.], length(data.tick.label.loc))
                #
                #horizontal dark grid lines  tick.label.loc
                #
                #
                matlines(rbind(pp.data(data.tick.label.loc, log.of.data), pp.data(data.tick.label.loc,
                        log.of.data)), rbind(yvec.low, yvec.high), col = 1., lty = 1., lwd = 1.)
                xvec.low <- rep(usr.out[1.], length(tick.label.loc))
                xvec.high <- rep(usr.out[2.], length(tick.label.loc))
                matlines(rbind(xvec.low, xvec.high), rbind(pp.quant(tick.label.loc, distribution, shape),
                        pp.quant(tick.label.loc, distribution, shape)), col = 1., lty = 1., lwd = 1.)
        }
        if(!slope.axis) {
                mtext(side = 2., line = 0, text = ylab, cex = cex.labs)
        }
        else {
                #
                #set up slope scale
                #
                add.slope.scale(distribution, shape, draw.line)
        }
        #
        #
        #we only need to adjust if in R
        # we will have different adjustments, depending on whether this is a 6-plot or not
        #
        use.title.line.adj <- title.line.adj
        #if(!checkIfR()) {
        #        #  here we are in S-PLUS
        #        use.title.line.adj <- 0
        #}
        #else {
                # here we are in R
                use.title.line.adj <- title.line.adj - 2
        #}
        #cat("*********************use.title.line.adj=", use.title.line.adj, "\n")
        #if(is.postscriptok() && substring(probplot.setup.title.out$new.title, 1., 1.) == "~") {
        #        mixed.mtext(side = 3., line = probplot.setup.title.out$title.line + use.title.line.adj, outer = F,
        #                texts = probplot.setup.title.out$new.title, adj = 0.5, cex = cex)
        #}
        #else {
        #        mtext(side = 3., line = probplot.setup.title.out$title.line + use.title.line.adj, outer = F,
        #                text = probplot.setup.title.out$new.title, cex = cex.title)
        #}
        #CheckPrintDataName()
        return(log.of.data)
}
################################################################################
GetAxesRange=function(type, x.axis, x.range, xlab, y.axis, y.range, ylab)
{
        #
        #if is.LockAxes(), try to recover
        #current axes on frame 0 for use later
        #
        #
        recovered <- F
        switch(type,
                probplot.setup = {
                        file.name <- ".axes.probplot.setup"
                }
                ,
                event.plot.setup = {
                        file.name <- ".axes.event.plot.setup"
                }
                ,
                plot.paper = {
                        file.name <- ".axes.plot.paper"
                }
                ,
                {
                        wqm.stop(paste(type, "is not recognized\n"))
                }
                )
        #
        #
        #
        #if  is.LockAxes(), the file is there, and there is a match, use it
        #otherwise send back the input
        #
       the.axes <- list(x.range = x.range, x.axis = x.axis, xlab = xlab, y.range = y.range, y.axis = y.axis,
                        ylab = ylab)
        #
        #
        #save whatever we have
        #


        return(the.axes)
}
################################################################################
get.prob.scales=function(distribution, ylab = GetSplidaDefault("SPLIDA.LabelOnYaxis"), prob.range = NULL, shape = NULL)
{
        sev.tick.location <- c(".00000001", ".00000002", ".00000003", ".00000005", ".0000001", ".0000002",
                ".0000003", ".0000005", ".000001", ".000002", ".000003", ".000005", ".00001", ".000012",
                ".000014", ".000016", ".000018", ".00002", ".000025", ".00003", ".000035", ".00004", ".000045",
                ".00005", ".00006", ".00007", ".00008", ".00009", ".0001", ".00012", ".00014", ".00016", ".00018",
                ".0002", ".00025", ".0003", ".0004", ".00045", ".0005", ".0006", ".0007", ".0008", ".0009",
                ".001", ".0012", ".0014", ".0016", ".0018", ".002", ".0025", ".003", ".0035", ".004", ".0045",
                ".005", ".006", ".007", ".008", ".009", ".01", ".012", ".014", ".016", ".018", ".02", ".022",
                ".024", ".026", ".028", ".03", ".035", ".04", ".045", ".05", ".06", ".07", ".08", ".09", ".1",
                ".12", ".14", ".16", ".18", ".2", ".22", ".24", ".26", ".28", ".3", ".35", ".4", ".45", ".5",
                ".55", ".6", ".65", ".7", ".75", ".8", ".85", ".9", ".92", ".94", ".96", ".98", ".99", ".995",
                ".999", ".9999", ".99999", ".999999", ".9999999")
        sev.tick.labels <- c(".00000001", ".00000002", ".00000003", ".00000005", ".0000001", ".0000002",
                ".0000003", ".0000005", ".000001", ".000002", ".000003", ".000005", ".00001", ".00002", ".00003",
                ".00005", ".0001", ".0002", ".0003", ".0005", ".001", ".003", ".005", ".01", ".02", ".03", ".05",
                ".1", ".2", ".3", ".5", ".7", ".9", ".98", ".99", ".999", ".9999", ".99999", ".999999",
                ".9999999")
        sev.Percent.tick.labels <- c(".000001", ".000002", ".000003", ".000005", ".00001", ".00002", ".00003",
                ".00005", ".0001", ".0002", ".0003", ".0005", ".001", ".002", ".003", ".005", ".01", ".02", ".03",
                ".05", ".1", ".3", ".5", "1", "2", "3", "5", "10", "20", "30", "50", "70", "90", "98", "99",
                "99.9", "99.99", "99.999", "99.9999", "99.99999")
        normal.tick.location <- c(".00000001", ".00000002", ".00000003", ".00000005", ".0000001", ".0000002",
                ".0000003", ".0000005", ".000001", ".000002", ".000003", ".000005", ".00001", ".000012",
                ".000014", ".000016", ".000018", ".00002", ".000025", ".00003", ".000035", ".00004", ".000045",
                ".00005", ".00006", ".00007", ".00008", ".00009", ".0001", ".0002", ".0003", ".0004", ".0005",
                ".0006", ".0007", ".0008", ".0009", ".001", ".0015", ".002", ".003", ".004", ".005", ".0075",
                ".01", ".015", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".1", ".12", ".14", ".16",
                ".18", ".2", ".25", ".3", ".35", ".4", ".45", ".5", ".55", ".6", ".65", ".7", ".75", ".8", ".82",
                ".84", ".86", ".88", ".9", ".91", ".92", ".93", ".94", ".95", ".96", ".97", ".98", ".985", ".99",
                ".9925", ".995", ".996", ".997", ".998", ".9985", ".999", ".9991", ".9992", ".9993", ".9994",
                ".9995", ".9996", ".9997", ".9998", ".9999", ".99999", ".999999", ".9999999")
        normal.tick.labels <- c(".00000001", ".00000002", ".00000003", ".00000005", ".0000001", ".0000002",
                ".0000003", ".0000005", ".000001", ".000002", ".000003", ".000005", ".00001", ".00005", ".0001",
                ".0005", ".001", ".002", ".005", ".01", ".02", ".05", ".1", ".2", ".3", ".4", ".5", ".6", ".7",
                ".8", ".9", ".95", ".98", ".99", ".995", ".998", ".999", ".9995", ".9999", ".99999", ".999999",
                ".9999999")
        normal.Percent.tick.labels <- c(".000001", ".000002", ".000003", ".000005", ".00001", ".00002", ".00003",
                ".00005", ".0001", ".0002", ".0003", ".0005", ".001", ".005", ".01", ".05", ".1", ".2", ".5",
                "1", "2", "5", "10", "20", "30", "40", "50", "60", "70", "80", "90", "95", "98", "99", "99.5",
                "99.8", "99.9", "99.95", "99.99", "99.999", "99.9999", "99.99999")
        exponential.tick.location <- c(".001", ".005", ".01", "0.05", ".1", ".15", ".2", ".25", ".30", ".35",
                ".4", ".45", ".4", ".5", ".55", ".6", ".45", ".65", ".7", ".75", ".8", ".85", ".9", ".91", ".92",
                ".93", ".94", ".95", ".96", ".97", ".98", ".985", ".991", ".992", ".993", ".994", ".995", ".996",
                ".997", ".998", ".9985", ".999", ".9991", ".9992", ".9993", ".9994", ".9995", ".9997", ".9998",
                ".9999", ".99999", ".9999999", ".9999999")
        exponential.tick.labels <- c(".001", ".01", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", ".95",
                ".98", ".99", ".995", ".998", ".999", ".9995", ".9998", ".9999", ".99999", ".999999", ".9999999")
        exponential.Percent.tick.labels <- c(".1", "1", "10", "20", "30", "40", "50", "60", "70", "80", "90",
                "95", "98", "99", "99.5", "99.8", "99.9", "99.95", "99.98", "99.99", "99.999", "99.9999",
                "99.99999")
        uniform.tick.location <- c("0", ".01", ".02", ".03", ".04", ".05", ".06", ".07", ".08", ".09", ".1",
                ".11", ".12", ".13", ".14", ".15", ".16", ".17", ".18", ".19", ".2", ".21", ".22", ".23", ".24",
                ".25", ".26", ".27", ".28", ".29", ".3", ".31", ".32", ".33", ".34", ".35", ".36", ".37", ".38",
                ".39", ".4", ".41", ".42", ".43", ".44", ".45", ".46", ".47", ".48", ".49", ".5", ".51", ".52",
                ".53", ".54", ".55", ".56", ".57", ".58", ".59", ".6", ".61", ".62", ".63", ".64", ".65", ".66",
                ".67", ".68", ".69", ".7", ".71", ".72", ".73", ".74", ".75", ".76", ".77", ".78", ".79", ".8",
                ".81", ".82", ".83", ".84", ".85", ".86", ".87", ".88", ".89", ".9", ".91", ".92", ".93", ".94",
                ".95", ".96", ".97", ".98", ".99", "1.")
        uniform.tick.labels <- c("0", ".1", ".2", ".3", ".4", ".5", ".6", ".7", ".8", ".9", "1.")
        uniform.Percent.tick.labels <- c("0", "10", "20", "30", "40", "50", "60", "70", "80", "90", "100")
        distribution <- generic.distribution(distribution)
        if(is.null(prob.range))
                prob.range <- range(as.numeric(sev.tick.labels))
        dist.info.out <- dist.info(distribution)
        #
        #check to make sure that we have the required shape parameters (if any)
        #
        if(length(shape) < dist.info.out$num.shape.needed) stop(paste("Need", dist.info.out$num.shape.needed,
                        "shape parameters for", distribution, "distribution--", length(shape), "provided"))
        #
        #get labels and locations according to distribbution
        #
        if(is.null(shape)) {
                switch(distribution,
                        gng = {
                                tick.labels <- sev.tick.labels
                                tick.location <- sev.tick.location
                                percent.tick.labels <- sev.Percent.tick.labels
                        }
                        ,
                        weibull = ,
                        sev = {
                                tick.labels <- sev.tick.labels
                                tick.location <- sev.tick.location
                                percent.tick.labels <- sev.Percent.tick.labels
                        }
                        ,
                        loglogistic = ,
                        logistic = {
                                tick.labels <- normal.tick.labels
                                tick.location <- normal.tick.location
                                percent.tick.labels <- normal.Percent.tick.labels
                        }
                        ,
                        lognormal = ,
                        normal = {
                                tick.labels <- normal.tick.labels
                                tick.location <- normal.tick.location
                                percent.tick.labels <- normal.Percent.tick.labels
                        }
                        ,
                        frechet = ,
                        lev = {
                                tick.labels <- normal.tick.labels
                                tick.location <- normal.tick.location
                                percent.tick.labels <- normal.Percent.tick.labels
                        }
                        ,
                        uniform = {
                                tick.labels <- wqm.pretty(prob.range, n = 5.)
                                tick.location <- wqm.pretty(tick.labels, n = 25.)
                                percent.tick.labels <- tick.labels * 100
                        }
                        ,
                        loguniform = {
                                tick.labels <- wqm.pretty(prob.range, n = 5.)
                                tick.location <- wqm.pretty(tick.labels, n = 25.)
                                percent.tick.labels <- wqm.pretty(prob.range, n = 5.) * 100
                        }
                        ,
                        exponential = {
                                if(prob.range[2.] > 0.01) {
                                        tick.labels <- exponential.tick.labels
                                        tick.location <- exponential.tick.location
                                        percent.tick.labels <- exponential.Percent.tick.labels
                                }
                                else {
                                        tick.labels <- sev.tick.labels
                                        tick.location <- sev.tick.location
                                        percent.tick.labels <- sev.Percent.tick.labels
                                }
                        }
                        ,
                        stop("Distribution not recognized"))
                if(dist.info.out$take.logs == "never")
                        logger <- F
                else logger <- T
        }
        else {
                #
                #here we have a specified shape parameter or parameters
                #
                if(!is.logdist(distribution)) stop(paste("Shape parameter provided for", distribution,
                                "distribution"))
                if(prob.range[2.] > 0.01) {
                        tick.location <- exponential.tick.location
                        tick.labels <- exponential.tick.labels
                        percent.tick.labels <- as.character(as.numeric(tick.labels) * 100)
                }
                else {
                        tick.location <- sev.tick.location
                        tick.labels <- sev.tick.labels
                        percent.tick.labels <- as.character(as.numeric(tick.labels) * 100)
                }
                if(dist.info.out$take.logs == "always")
                        logger <- T
                else logger <- F
        }
        return(list(tick.location = tick.location, tick.labels = tick.labels, percent.tick.labels =
                percent.tick.labels, logger = logger, distribution = dist.info.out$formal.name, prob.scale =
                dist.info.out$prob.scale))
}
################################################################################
generic.distribution=function(distribution, allow = F)
{
        #
        #if allow=T, the this can be used as a check to see if a char string
        #       is a valid distribution--returns NULL in this case
        #       otherwise there is a stop
        #
        #
        if(!is.character(distribution)) stop(paste("distribution must be character string:", distribution))
        switch(distribution,
                Uniform = ,
                uniform = {
                        distribution <- "uniform"
                }
                ,
                "Log-Uniform" = ,
                LogUniform = ,
                loguniform = {
                        distribution <- "loguniform"
                }
                ,
                exponential = ,
                Exponential = {
                        distribution <- "exponential"
                }
                ,
                "Smallest Extreme Value" = ,
                "Smallest extreme value" = ,
                Smallestextremevalue = ,
                sev = ,
                Sev = {
                        distribution <- "sev"
                }
                ,
                logsev = ,
                LogSev = ,
                Logsev = ,
                Weibull = ,
                weibull = {
                        distribution <- "weibull"
                }
                ,
                nor = ,
                Normal = ,
                normal = {
                        distribution <- "normal"
                }
                ,
                "Log-Normal" = ,
                "Log-normal" = ,
                LogNormal = ,
                Lognormal = ,
                lognormal = ,
                Lognormal.basee = ,
                lognormal.basee = {
                        distribution <- "lognormal"
                }
                ,
                Lognormal10 = ,
                Lognormal.base10 = ,
                lognormal.base10 = {
                        distribution <- "lognormal10"
                }
                ,
                LogLogistic = ,
                Loglogistic = ,
                loglogistic = {
                        distribution <- "loglogistic"
                }
                ,
                Logistic = ,
                logistic = {
                        distribution <- "logistic"
                }
                ,
                "Largest Extreme Value" = ,
                "Largest extreme value" = ,
                Largestextremevalue = ,
                lev = {
                        distribution <- "lev"
                }
                ,
                gamma = ,
                Gamma = {
                        distribution <- "gamma"
                }
                ,
                IGAU = ,
                Igau = ,
                igau = ,
                "Inverse Gaussian" = {
                        distribution <- "igau"
                }
                ,
                BISA = ,
                Bisa = ,
                bisa = ,
                "Birnbaum-Saunders" = {
                        distribution <- "bisa"
                }
                ,
                GOMA = ,
                Goma = ,
                goma = ,
                "Gompertz-Makeham" = {
                        distribution <- "goma"
                }
                ,
                "Generalized Gamma" = ,
                GNG = ,
                gng = {
                        distribution <- "gng"
                }
                ,
                "Generalized F" = ,
                gnf = {
                        distribution <- "gnf"
                }
                ,
                "Log Extended Generalized Gamma" = ,
                egengl = {
                        distribution <- "egengl"
                }
                ,
                "Extended Generalized Gamma" = ,
                egeng = {
                        distribution <- "egeng"
                }
                ,
                "Reciprocal Weibull" = ,
                Frechet = ,
                frechet = ,
                "Log Lev" = ,
                loglev = ,
                "recip_weibull" = ,
                recipweibull = {
                        distribution <- "frechet"
                }
                ,
                "Sev Generalized Threshold Scale" = ,
                sevgets = {
                        distribution <- "sevgets"
                }
                ,
                "Lev Generalized Threshold Scale" = ,
                levgets = {
                        distribution <- "levgets"
                }
                ,
                "Normal Generalized Threshold Scale" = ,
                norgets = ,
                normalgets = {
                        distribution <- "normalgets"
                }
                ,
                Beta = ,
                beta = {
                        distribution <- "beta"
                }
                ,
                "Log-beta" = ,
                Logbeta = ,
                "Log-Beta" = ,
                logbeta = {
                        distribution <- "logbeta"
                }
                ,
                Triangle = ,
                triangle = {
                        distribution <- "triangle"
                }
                ,
                "Log-Triangle" = ,
                "Log-triangle" = ,
                Logtriangle = ,
                logtriangle = {
                        distribution <- "logtriangle"
                }
                ,
                xdummy = {
                        distribution <- "xdummy"
                }
                ,
                {
                        #
                        #
                        #browser()
                        #null return needed for some checking in summary
                        #
                        #
                        if(allow) {
                                return(NULL)
                        }
                        else {
                                stop(paste(distribution,
                                        "is an unrecognized distribution in generic.distribution()"))
                        }
                }
                )
        return(distribution)
}
################################################################################
dist.info=function(distribution, allow = F)
{
        #
        #distribution could be null if it is not in the list
        #
        collapse.distribution <- paste(distribution, collapse = ",")
        distribution <- generic.distribution(collapse.distribution, allow = allow)
        if(is.null(distribution)) {
                #
                #for special, unnumbered distributions
                #
                take.logs <- "never"
                num.shape.needed <- 0.
                formal.name <- collapse.distribution
                shape <- F
                prob.scale <- "sev"
                idist <- 0
                return(idist, take.logs, num.shape.needed, formal.name)
        }
        idist <- numdist(distribution)
        switch(distribution,
                sev = {
                        take.logs <- "never"
                        num.shape.needed <- 0.
                        formal.name <- "Smallest Extreme Value"
                        shape <- F
                        prob.scale <- "sev"
                }
                ,
                weibull = {
                        take.logs <- "if.no.shape"
                        num.shape.needed <- 0.
                        formal.name <- "Weibull"
                        prob.scale <- "sev"
                }
                ,
                uniform = {
                        take.logs <- "never"
                        num.shape.needed <- 0.
                        formal.name <- "Uniform"
                        prob.scale <- "uniform"
                }
                ,
                loguniform = {
                        take.logs <- "always"
                        num.shape.needed <- 0.
                        formal.name <- "Log-Uniform"
                        prob.scale <- "loguniform"
                }
                ,
                normal = {
                        take.logs <- "never"
                        num.shape.needed <- 0.
                        formal.name <- "Normal"
                        prob.scale <- "normal"
                }
                ,
                lognormal = {
                        take.logs <- "if.no.shape"
                        num.shape.needed <- 0.
                        formal.name <- "Lognormal"
                        prob.scale <- "normal"
                }
                ,
                logistic = {
                        take.logs <- "never"
                        num.shape.needed <- 0.
                        formal.name <- "Logistic"
                        prob.scale <- "logistic"
                }
                ,
                loglogistic = {
                        take.logs <- "if.no.shape"
                        num.shape.needed <- 0.
                        formal.name <- "Loglogistic"
                        prob.scale <- "logistic"
                }
                ,
                exponential = {
                        take.logs <- "never"
                        num.shape.needed <- 0.
                        formal.name <- "Exponential"
                        prob.scale <- "exponential"
                }
                ,
                gamma = {
                        take.logs <- "never"
                        num.shape.needed <- 1.
                        formal.name <- "Gamma"
                        prob.scale <- "gamma"
                }
                ,
                gng = {
                        take.logs <- "always"
                        num.shape.needed <- 1.
                        formal.name <- "Generalized Gamma"
                        prob.scale <- "gen-gamma"
                }
                ,
                lev = {
                        take.logs <- "never"
                        num.shape.needed <- 0.
                        formal.name <- "Largest Extreme Value"
                        prob.scale <- "normal"
                }
                ,
                frechet = {
                        take.logs <- "if.no.shape"
                        num.shape.needed <- 0.
                        formal.name <- "Frechet"
                        prob.scale <- "normal"
                }
                ,
                igau = {
                        take.logs <- "never"
                        num.shape.needed <- 1.
                        formal.name <- "Inverse Gaussian"
                        prob.scale <- "igau"
                }
                ,
                bisa = {
                        take.logs <- "never"
                        num.shape.needed <- 1.
                        formal.name <- "Birnbaum-Saunders"
                        prob.scale <- "bisa"
                }
                ,
                goma = {
                        take.logs <- "never"
                        num.shape.needed <- 2.
                        formal.name <- "Gompertz-Makeham"
                        prob.scale <- "goma"
                }
                ,
                gnf = {
                        take.logs <- "always"
                        num.shape.needed <- 2.
                        formal.name <- "Generalized F"
                        prob.scale <- "gen-F"
                }
                ,
                normalgets = {
                        take.logs <- "never"
                        num.shape.needed <- 1.
                        formal.name <- "Normal Generalized Threshold Scale"
                        prob.scale <- "lognormal"
                }
                ,
                sevgets = {
                        take.logs <- "never"
                        num.shape.needed <- 1.
                        formal.name <- "Smallest Extreme Value Generalized Threshold Scale"
                        prob.scale <- "weibull"
                }
                ,
                egeng = {
                        take.logs <- "always"
                        num.shape.needed <- 1.
                        formal.name <- "Extended Generalized Gamma"
                        prob.scale <- "weibull"
                }
                ,
                {
                        if(allow)
                                return(list())
                        stop("Distribution not recognized")
                }
                )
        return(list(idist = idist, take.logs = take.logs, num.shape.needed = num.shape.needed, formal.name =
                formal.name))
}

################################################################################
numdist=function(distribution, allow = F)
{
        #
        #return code number for a distribution
        #
        the.distribution <- generic.distribution(distribution, allow = allow)
        if(is.null(the.distribution) || distribution == "") {
                the.message <- paste("Distribution not recognized in numdist:", distribution)
                if(allow) {
                        idist <- 0.
                }
                else {
                        stop(the.message)
                }
        }
        else {
                switch(generic.distribution(distribution, allow = allow),
                        sev = idist <- 1.,
                        weibull = idist <- 2.,
                        normal = idist <- 3.,
                        lognormal = idist <- 4.,
                        logistic = idist <- 5.,
                        loglogistic = idist <- 6.,
                        lev = idist <- 7.,
                        frechet = idist <- 8.,
                        gng = idist <- 10.,
                        gamma = idist <- 12.,
                        logexponential = idist <- 13.,
                        exponential = idist <- 14.,
                        igau = idist <- 16.,
                        bisa = idist <- 18.,
                        goma = idist <- 20.,
                        gnf = idist <- 22.,
                        uniform = idist <- 23.,
                        loguniform = idist <- 24.,
                        egeng = idist <- 26.,
                        sevgets = idist <- 27.,
                        normalgets = idist <- 29.,
                        levgets = idist <- 30.,
                        beta = idist <- 31.,
                        logbeta = idist <- 32.,
                        triangle = idist <- 33.,
                        logtriangle = idist <- 34.,
                        )
        }
        return(idist)
}

################################################################################
pp.quant=function(p, distribution, shape = NULL)
{
        the.attributes <- attributes(p)
        if(is.null(shape) || distribution == "uniform")
                answer <- quant(p, distribution)
        else answer <- gquant(p, distribution, shape)
        attributes(answer) <- the.attributes
        return(answer)
}

################################################################################
quant=function(p, distribution)
{
        #
        #this version of quant assumes that we are transforming to a log scale
        #
        #see gquant for general quant function, by name and abbrev
        #
        if(any(p<0)  | any(p>1) ) error("quantiles must be between 0 and 1")
        switch(generic.distribution(distribution),
                weibull = ,
                sev = {
                        qvec <- logb( - logb(1. - p))
                }
                ,
                uniform = ,
                loguniform = {
                        qvec <- p
                }
                ,
                frechet = ,
                lev = {
                        qvec <-  - logb( - logb(p))
                }
                ,
                normal = ,
                lognormal = {
                        qvec <- qnorm(p)
                }
                ,
                loglogistic = ,
                logistic = {
                        qvec <- qlogis(p)
                }
                ,
                exponential = {
                        qvec <- qexp(p)
                }
                ,
                stop(paste(distribution, "is unrecognized distribution in quant()")))
        return(qvec)
}

################################################################################
logax=function(xmin, xmax, nint = 5., ntick = 4., which.labels = NULL, ...)
{
        if(missing(xmax) && length(xmin) == 2.) {
                xmax <- xmin[2.]
                xmin <- xmin[1.]
        }
        if(is.null(xmax)) {
                wqm.warning("xmax NULL in logax; set xmax <- .1")
                xmax <- 0.10000000000000001
        }
        if(log10(xmax) >= Inf) {
                wqm.warning("xmax is too big, set equal to 10^300")
                xmax <- 10.^300.
        }
        if(xmax <= 0.)
                stop(paste("non-positive xmax=", xmax, "in logax"))
        if(xmax < xmin)
                wqm.warning(paste("log axis xmax=", xmax, "< xmin=", xmin, "in logax"))
        if(xmin <= 0.) {
                xmin.fix <- min(xmax/100., 1.0000000000000001e-005)
                wqm.warning(paste("negative xmin=", xmin, "with xmax=", xmax, "min set to", xmin.fix))
                xmin <- xmin.fix
        }
        drange <- log10(xmax) - log10(xmin)
        if(drange > 300.) {
                drange <- log10(xmax) - log10(xmin)
                wqm.warning(paste("lxmax=", xmax, "xmin=", xmin, "dynamic range=", drange,
                        "is too large; restricted to 300"))
                xmean <- ((log10(xmax) + log10(xmin))/2.)
                xmax <- 10.^(xmean + 150.)
                xmin <- 10.^(xmean - 150.)
                cat(paste("\nSet to xmax=", xmax, "xmin=", xmin, "\n"))
        }
        drange <- xmax/xmin
        num.dec <- floor(log10(drange)) + 4.
        minmult <- floor(log10(xmin)) - 1.
        mult.vec <- 10.^seq(minmult, minmult + num.dec - 1.)
        if(is.null(which.labels))
                which.labels <- num.dec - 1.
        else {
                switch(which.labels,
                        which.labels <- 1.,
                        which.labels <- 5.,
                        which.labels <- 10.)
        }
        switch(as.character(which.labels),
                "0" = ,
                "1" = ,
                "2" = ,
                "3" = {
                        linax.out <- linax(xmax, xmin, nint = 5., ntick = ntick)
                        numlab <- as.numeric(linax.out$ticlab)
                        if(any(numlab > 0.)) {
                                linax.out <- linax(xmax, xmin, nint = 10., ntick = ntick)
                                numlab <- as.numeric(linax.out$ticlab)
                        }
                        numloc <- as.numeric(linax.out$ticloc)
                        return(list(ticlab = format(linax.out$ticlab[numlab > 0.]), ticloc = format(linax.out$
                                ticloc[numloc > 0.])))
                }
                ,
                "5" = ,
                "4" = {
                        ticlab <- rep(c(1., 2., 5.), num.dec) * sort(rep(mult.vec, 3.))
                        ticloc <- rep(c(1., 1.2, 1.3999999999999999, 1.6000000000000001, 1.8, 2., 3., 4., 5.,
                                6., 7., 8., 9.), num.dec) * sort(rep(mult.vec, 13.))
                }
                ,
                {
                        ticlab <- rep(c(1.), num.dec) * sort(rep(mult.vec, 1.))
                        ticloc <- rep(c(1., 2., 3., 4., 5., 6., 7., 8., 9.), num.dec) * sort(rep(mult.vec, 9.))
                }
                )
        xminp <- max(ticlab[ticlab <= xmin])
        xmaxp <- min(ticlab[ticlab >= xmax])
        in.range.lab <- ticlab >= xminp & ticlab <= xmaxp
        ticlab <- ticlab[in.range.lab]
        in.range.loc <- ticloc >= min(ticlab) & ticloc <= max(ticlab)
        ticloc <- ticloc[in.range.loc]
        return(list(ticlab = format(ticlab), ticloc = format(ticloc)))
}

################################################################################
linax=function(xmax, xmin, nint = 5., nticks = GetSplidaDefault("SPLIDA.NumberTicks"))
{
        #
        #       major and minor tic label locations
        #
        #       nint controls number of major label locations
        #       nticks controls number of ticks within label locations
        #
        if(missing(xmin) && length(xmax == 2.)) {
                xmin <- xmax[1.]
                xmax <- xmax[2.]
        }
        ticlab <- wqm.pretty(c(xmax, xmin), n = nint)
        ticloc <- seq(ticlab[1.], ticlab[length(ticlab)], length = (length(ticlab) - 1.) * (nticks + 1.) + 1.)
        return(list(ticlab = format(ticlab), ticloc = format(ticloc)))
}

################################################################################
wqm.pretty=function(x, n = 5, min.n = n %/% 3, shrink.sml = 0.75, high.u.bias = 1.5, u5.bias = 0.5 + 1.5 * high.u.bias,
        eps.correct = 0, nint = 5)
{
        #
        #This is a version of pretty that should work in either S-PLUS or R
        #
        sdfd=1
        if(sdfd==1) {
                #here we are in R
                if(missing(nint)) pretty(x = x, n = n, min.n = min.n, shrink.sml = shrink.sml, high.u.bias = 1.5,
                                u5.bias = high.u.bias, eps.correct = eps.correct) else {
                        pretty(x = x, n = nint, min.n = min.n, shrink.sml = shrink.sml, high.u.bias = 1.5,
                                u5.bias = high.u.bias, eps.correct = eps.correct)
                }
        }
        else {
                pretty(x = x, nint = nint)
        }
}

################################################################################
pp.data=function(data.vector, log.of.data)
{
        if(log.of.data)
                return(logb(data.vector))
        else return(data.vector)
}

################################################################################
fix.exp.labels=function(ticlabels)
{
        ticlabels.orig <- ticlabels
        nchar.out <- nchar(ticlabels)
        where.zero <- rep(0., length(nchar.out))
        plus.or.neg <- rep("0", length(nchar.out))
        mantissa <- rep(NA, length(nchar.out))
        e.here <- rep(NA, length(nchar.out))
        power <- rep(NA, length(nchar.out))
        charmat <- matrix(NA, nrow = length(nchar.out), ncol = max(nchar.out))
        ticlabels <- paste(ticlabels, "     ")
        for(i in 1.:max(nchar.out)) {
                #check to see if we have hit them all yet
                charmat[, i] <- substring(ticlabels, i, i)
                e.here[charmat[, i] == "e"] <- i
        }
        strip.lead0 <- function(xstring)
        {
                nchar.out <- nchar(xstring)
                fchar <- substring(xstring, 1., 1.)
                xstring <- ifelse(fchar == "0", substring(xstring, 2., nchar.out), xstring)
                return(xstring)
        }
        worknow <- !is.na(e.here)
        if(any(worknow)) {
                plus.or.neg[worknow] <- substring(ticlabels[worknow], e.here[worknow] + 1., e.here[worknow] +
                        1.)
                mantissa[worknow] <- substring(ticlabels[worknow], 1., e.here[worknow] - 1.)
                power[worknow] <- strip.lead0(substring(ticlabels[worknow], e.here[worknow] + 2., nchar.out[
                        worknow]))
                was.exp <- !is.na(mantissa)
                man1 <- mantissa == "1"
                mantissa[man1] <- ""
                mantissa[!man1] <- paste(mantissa[!man1], "x", sep = "")
        }
        was.exp <- !is.na(mantissa)
        #
        #need to fix this up with mixed.mtext to replace axes, if possible
        #if(names(dev.cur()) == "postscript"){
        #ticlabels.orig[was.exp ] <-
        #  paste("~.",mantissa[was.exp ],"10~u.8~c.7~.",plus.or.neg[was.exp],power[was.exp ],"~",sep="")}
        #else{
        plus.or.neg[plus.or.neg != "-"] <- ""
        ticlabels.orig[was.exp] <- paste(mantissa[was.exp], "10", "^", plus.or.neg[was.exp], power[was.exp],
                sep = "")
        return(ticlabels.orig)
}

################################################################################
vector.power10=function(str.vec)
{
        retvec <- rep(NA, length(str.vec))
        for(i in seq(along = str.vec)) {
                retvec[i] <- power10(str.vec[i])
        }
        return(retvec)
}

################################################################################
power10=function(xstring, maxlen = 10.)
{
        parse.vector <- substring(substring(xstring, 1.:maxlen), 1., 1.)
        e.pos <- seq(1., maxlen)[parse.vector == "e"]
        if(length(e.pos) == 0.) {
                #strip blanks and nulls
                #remove leading 0
                xstring <- strip.blanks.nulls(xstring)
                if(substring(xstring, 1., 1.) == 0.)
                        xstring <- substring(xstring, 2.)
                return(xstring)
        }
        null.pos <- min(seq(1., maxlen)[parse.vector == ""])
        mantissa <- substring(xstring, 1., e.pos - 1.)
        #remove leading + from power
        power <- substring(xstring, e.pos + 1., null.pos - 1.)
        #remove leading 0 from power
        if(substring(power, 1., 1.) == "+") power <- substring(power, 2.)
        #remove leading -0 from power
        if(substring(power, 1., 1.) == "0") power <- substring(power, 2.)
        if(substring(power, 2., 2.) == "0" && substring(power, 1., 1.) == "-")
                power <- paste("-", substring(power, 3.), sep = "", collapse = "")
        if(mantissa == "1") {
                mantissa <- ""
                xnumber <- paste("~.10~u.8~c.8~.", power, "~c1~.", sep = "")
        }
        else {
                xnumber <- paste("~.", mantissa, "x10~u.8~c.8~.", power, "~c1~.", sep = "")
        }
        return(xnumber)
}

################################################################################
strip.blanks.nulls=function(xstring, maxlen = nchar(xstring))
{
        parse.vector <- substring(substring(xstring, 1.:maxlen), 1., 1.)
        blanks <- parse.vector == " " | parse.vector == ""
        return(paste(parse.vector[!blanks], collapse = ""))
}

################################################################################
blank.fix=function(ss)
{
  nn=length(ss)
  for (i in 1:nn)
  {
    tt=ss[i]
    mm=nchar(tt)
    ww=NULL
    for(j in 1:mm)
    {
      ee=substr(tt,j,j)
      if(ee!=" ")
      {
        ww=paste(ww,ee,sep="")
      }
    }
    ss[i]=ww
  }
  return(ss)
}

################################################################################
multi.norm.sim=function(n,mu,sigma,scale=1)
{
  #save(sigma,file="tmp.sigma")
  p=dim(sigma)[1]
  res=matrix(rnorm(n*p),nrow=p,ncol=n)
  res=res*scale

  aa=eigen(sigma,symmetric=T)

  #uu=diag(aa$d,length(aa$d),length(aa$d))
  #dd=aa$v%*%sqrt(uu)

  #print(dim(aa$vectors))
  #print(dim(diag(sqrt(aa$values),)))

  dd=aa$vectors%*%diag(sqrt(aa$values),)
  
  
  res=dd%*%res

  res=sweep(res,1,mu,"+")
  
  res=t(res)
  
  return(res)
}
################################################################################
m.spline.x=function(x,tt,i,k)
{
  ti=tt[i]
  tik=tt[i+k]
  if(x<ti|x>=tik)#x<ti | x>tik |ti==tik)
  {
     res=0
     return(res)
  }else{
        if(k==1)
        {
          res=1/(tik-ti)
          return(res)
        }else{
              a1=m.spline.x(x,tt,i,k-1)
              a2=m.spline.x(x,tt,i+1,k-1)
              res=(k*((x-ti)*a1+(tik-x)*a2))/((k-1)*(tik-ti))
              return(res)
             }
       }
}
################################################################################
i.spline.x=function(xx,tt,i,k,delta=0.0001,Cs=F)
{
  nn=length(xx)
  a1=min(tt)
  b1=max(tt)
  xtt=seq(a1,b1,,by=delta)
  nxtt=length(xtt)
  ytt=double(nxtt)
  for(j in 1:nxtt)
  {

    ytt[j]=m.spline.x(xtt[j],tt,i,k)
  }
  ytt=cumsum(ytt)*delta
  if(Cs)
  {
    ytt=cumsum(ytt)*delta
  }

  ii=ceiling((xx-a1)/delta)
  ii[ii>nxtt]=nxtt
  ii[ii==0]=1
  res=ytt[ii]
  return(res)
}
################################################################################
MIC.splines.basis=function(x, df = NULL, knots = NULL,
boundary.knots=NULL,type="Ms",degree = 3,delta=0.01,eq.alloc=F)
{
  if(is.null(df)&is.null(knots))
  {
    stop("either df or knots needs to be supplied")
  }
  
  if(!is.null(df))
  {
    cc=df-degree+1
    knots=quantile(x,1:(cc-1)/cc)
  }
  
  if(is.null(boundary.knots))
  {
   boundary.knots=range(x)
   boundary.knots[2]=boundary.knots[2]+0.0001
  }
  
  if(eq.alloc==T)
  {
    if(is.null(df))
    {
      stop("df needs to be supplied for equal allocation")
    }
    cc=df-degree+1
    knots=boundary.knots[1]+(1:(cc-1))*(boundary.knots[2]-boundary.knots[1])/cc
  }
  
  tt=c(rep(boundary.knots[1],degree),knots,rep(boundary.knots[2],degree))
  nn=length(x)
  mm=degree+length(knots)
  mat=matrix(0, nrow=nn,ncol=mm)
  
  #M splines basis
  if(type=="Ms")
  {
   for(i in 1:nn)
   {
    for(j in 1:mm)
    {
      mat[i,j]=m.spline.x(x[i],tt,j,k=degree)
    }
   }
  }

  #I spline
  if(type=="Is")
  {
   for(j in 1:mm)
   {
     mat[,j]=i.spline.x(x,tt,j,k=degree,delta=delta,Cs=F)
   }
   mat=sweep(mat,2,colMeans(mat),"-")
  }
  
  #C spines
  if(type=="Cs")
  {
   for(j in 1:mm)
   {
     mat[,j]=i.spline.x(x,tt,j,k=degree,delta=delta,Cs=T)
   }
   x1=x-mean(x)
   x1=x1/sqrt(sum(x1*x1))
   mat=sweep(mat,2,colMeans(mat),"-")
   for(j in 1:mm)
   {
     mat[,j]=mat[,j]-x1*sum(mat[,j]*x1)
   }
   mat=cbind(x1,mat)
  }
  ###return results
  res=list(mat=mat,x=x, df=df, knots=knots,boundary.knots=boundary.knots,
  type=type,degree=degree,delta=delta)
  attr(res,"class")="MICsplines"
  return(res)
}
################################################################################
plot.MICsplines=function(obj)
{
  mat=obj$mat
  x=obj$x
  matplot(x[order(x)],mat[order(x),],type="l",xlab="",ylab="",main=paste(obj$type,"Splines Basis"))
  
}
################################################################################
newlegend=function (x, y = NULL, legend, fill = NULL, col = par("col"),
    border = "black", lty, lwd, pch, angle = 45, density = NULL,
    bty = "o", bg = par("bg"), box.lwd = par("lwd"), box.lty = par("lty"),
    box.col = par("fg"), pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd,
    xjust = 0, yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0,
        0.5), text.width = NULL, text.col = par("col"), merge = do.lines &&
        has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE,
    title = NULL, inset = 0, xpd, title.col = text.col,aa=2)
{
    if (missing(legend) && !missing(y) && (is.character(y) ||
        is.expression(y))) {
        legend <- y
        y <- NULL
    }
    mfill <- !missing(fill) || !missing(density)
    if (!missing(xpd)) {
        op <- par("xpd")
        on.exit(par(xpd = op))
        par(xpd = xpd)
    }
    title <- as.graphicsAnnot(title)
    if (length(title) > 1)
        stop("invalid title")
    legend <- as.graphicsAnnot(legend)
    n.leg <- if (is.call(legend))
        1
    else length(legend)
    if (n.leg == 0)
        stop("'legend' is of length 0")
    auto <- if (is.character(x))
        match.arg(x, c("bottomright", "bottom", "bottomleft",
            "left", "topleft", "top", "topright", "right", "center"))
    else NA
    if (is.na(auto)) {
        xy <- xy.coords(x, y)
        x <- xy$x
        y <- xy$y
        nx <- length(x)
        if (nx < 1 || nx > 2)
            stop("invalid coordinate lengths")
    }
    else nx <- 0
    xlog <- par("xlog")
    ylog <- par("ylog")
    rect2 <- function(left, top, dx, dy, density = NULL, angle,
        ...) {
        r <- left + dx
        if (xlog) {
            left <- 10^left
            r <- 10^r
        }
        b <- top - dy
        if (ylog) {
            top <- 10^top
            b <- 10^b
        }
        rect(left, top, r, b, angle = angle, density = density,
            ...)
    }
    segments2 <- function(x1, y1, dx, dy, ...) {
        x2 <- x1 + dx
        if (xlog) {
            x1 <- 10^x1
            x2 <- 10^x2
        }
        y2 <- y1 + dy
        if (ylog) {
            y1 <- 10^y1
            y2 <- 10^y2
        }
        segments(x1, y1, x2, y2, ...)
    }
    points2 <- function(x, y, ...) {
        if (xlog)
            x <- 10^x
        if (ylog)
            y <- 10^y
        points(x, y, ...)
    }
    text2 <- function(x, y, ...) {
        if (xlog)
            x <- 10^x
        if (ylog)
            y <- 10^y
        text(x, y, ...)
    }
    if (trace)
        catn <- function(...) do.call("cat", c(lapply(list(...),
            formatC), list("\n")))
    cin <- par("cin")
    Cex <- cex * par("cex")
    if (is.null(text.width))
        text.width <- max(abs(strwidth(legend, units = "user",
            cex = cex)))
    else if (!is.numeric(text.width) || text.width < 0)
        stop("'text.width' must be numeric, >= 0")
    xc <- Cex * xinch(cin[1L], warn.log = FALSE)
    yc <- Cex * yinch(cin[2L], warn.log = FALSE)
    if (xc < 0)
        text.width <- -text.width
    xchar <- xc
    xextra <- 0
    yextra <- yc * (y.intersp - 1)
    ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
    ychar <- yextra + ymax
    if (trace)
        catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra,
            ychar))
    if (mfill) {
        xbox <- xc * 0.8
        ybox <- yc * 0.5
        dx.fill <- xbox
    }
    do.lines <- (!missing(lty) && (is.character(lty) || any(lty >
        0))) || !missing(lwd)
    n.legpercol <- if (horiz) {
        if (ncol != 1)
            warning("horizontal specification overrides: Number of columns := ",
                n.leg)
        ncol <- n.leg
        1
    }
    else ceiling(n.leg/ncol)
    has.pch <- !missing(pch) && length(pch) > 0
    if (do.lines) {
        x.off <- if (merge)
            -0.7
        else 0
    }
    else if (merge)
        warning("'merge = TRUE' has no effect when no line segments are drawn")
    if (has.pch) {
        if (is.character(pch) && !is.na(pch[1L]) && nchar(pch[1L],
            type = "c") > 1) {
            if (length(pch) > 1)
                warning("not using pch[2..] since pch[1L] has multiple chars")
            np <- nchar(pch[1L], type = "c")
            pch <- substr(rep.int(pch[1L], np), 1L:np, 1L:np)
        }
    }
    if (is.na(auto)) {
        if (xlog)
            x <- log10(x)
        if (ylog)
            y <- log10(y)
    }
    if (nx == 2) {
        x <- sort(x)
        y <- sort(y)
        left <- x[1L]
        top <- y[2L]
        w <- diff(x)
        h <- diff(y)
        w0 <- w/ncol
        x <- mean(x)
        y <- mean(y)
        if (missing(xjust))
            xjust <- 0.5
        if (missing(yjust))
            yjust <- 0.5
    }
    else {
        h <- (n.legpercol + (!is.null(title))) * ychar + yc
        w0 <- text.width + (x.intersp + 1) * xchar
        if (mfill)
            w0 <- w0 + dx.fill
        if (do.lines)
            w0 <- w0 + (aa + x.off) * xchar
        w <- ncol * w0 + 0.5 * xchar
        if (!is.null(title) && (abs(tw <- strwidth(title, units = "user",
            cex = cex) + 0.5 * xchar)) > abs(w)) {
            xextra <- (tw - w)/2
            w <- tw
        }
        if (is.na(auto)) {
            left <- x - xjust * w
            top <- y + (1 - yjust) * h
        }
        else {
            usr <- par("usr")
            inset <- rep(inset, length.out = 2)
            insetx <- inset[1L] * (usr[2L] - usr[1L])
            left <- switch(auto, bottomright = , topright = ,
                right = usr[2L] - w - insetx, bottomleft = ,
                left = , topleft = usr[1L] + insetx, bottom = ,
                top = , center = (usr[1L] + usr[2L] - w)/2)
            insety <- inset[2L] * (usr[4L] - usr[3L])
            top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3L] +
                h + insety, topleft = , top = , topright = usr[4L] -
                insety, left = , right = , center = (usr[3L] +
                usr[4L] + h)/2)
        }
    }
    if (plot && bty != "n") {
        if (trace)
            catn("  rect2(", left, ",", top, ", w=", w, ", h=",
                h, ", ...)", sep = "")
        rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
            lwd = box.lwd, lty = box.lty, border = box.col)
    }
    xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1),
        rep.int(n.legpercol, ncol)))[1L:n.leg]
    yt <- top - 0.5 * yextra - ymax - (rep.int(1L:n.legpercol,
        ncol)[1L:n.leg] - 1 + (!is.null(title))) * ychar
    if (mfill) {
        if (plot) {
            fill <- rep(fill, length.out = n.leg)
            rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox,
                col = fill, density = density, angle = angle,
                border = border)
        }
        xt <- xt + dx.fill
    }
    if (plot && (has.pch || do.lines))
        col <- rep(col, length.out = n.leg)
    if (missing(lwd))
        lwd <- par("lwd")
    if (do.lines) {
        seg.len <- aa
        if (missing(lty))
            lty <- 1
        lty <- rep(lty, length.out = n.leg)
        lwd <- rep(lwd, length.out = n.leg)
        ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
        if (trace)
            catn("  segments2(", xt[ok.l] + x.off * xchar, ",",
                yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
        if (plot)
            segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len *
                xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l],
                col = col[ok.l])
        xt <- xt + (seg.len + x.off) * xchar
    }
    if (has.pch) {
        pch <- rep(pch, length.out = n.leg)
        pt.bg <- rep(pt.bg, length.out = n.leg)
        pt.cex <- rep(pt.cex, length.out = n.leg)
        pt.lwd <- rep(pt.lwd, length.out = n.leg)
        ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
        x1 <- (if (merge && do.lines)
            xt - (seg.len/2) * xchar
        else xt)[ok]
        y1 <- yt[ok]
        if (trace)
            catn("  points2(", x1, ",", y1, ", pch=", pch[ok],
                ", ...)")
        if (plot)
            points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok],
                bg = pt.bg[ok], lwd = pt.lwd[ok])
    }
    xt <- xt + x.intersp * xchar
    if (plot) {
        if (!is.null(title))
            text2(left + w/2, top - ymax, labels = title, adj = c(0.5,
                0), cex = cex, col = title.col)
        text2(xt, yt, labels = legend, adj = adj, cex = cex,
            col = text.col)
    }
    invisible(list(rect = list(w = w, h = h, left = left, top = top),
        text = list(x = xt, y = yt)))
}


################################################################################
obj.size=function()
{
   a=ls(envir=.GlobalEnv)
   nn=length(a)
   ss=double(nn)
   for(i in 1:nn)
   {
     ss[i]=as.numeric(object.size(get(a[i])))/1024^2
   }
   res=data.frame(objs=a,size=ss)
   res=res[order(res[,2],decreasing=T),]
   return(res)
}

################################################################################
################################################################################