# ===================================================================
# === Ancilary indicators ===========================================
# ===================================================================

# Linear model version
slp2<-function(x,mat,ind){

  lm(y~x1,data.frame(x1=1:length(ind),y=log(mat[x,ind])))$coef[2]

}

# MLE version
slp<-function(x,mat,ind){

  x1<-1:length(ind)
  y=log(mat[x,ind])
  mux<-mean(x1,na.rm=T)
  muy<-mean(y,na.rm=T)
  SS<-sum((x1-mux)^2,na.rm=T)
  (1/SS)*sum((x1-mux)*(y-muy),na.rm=T)

}

# Average annual variability
AAV<-function(x,mat,ind){
  ni<-length(ind)
  mean(abs((mat[x,ind[2:ni]]-mat[x,ind[1:(ni-1)]])/mat[x,ind[1:(ni-1)]]))
}

# Log mean
mu<-function(x,mat,ind){
  log(mean(mat[x,ind]))
}

#' Characterize posterior predictive data
#'
#' @param PPD An object of class Data stored in the Misc slot of an MSE object following a call of \code{runMSE(PPD = TRUE)}.
#' @param styr Positive integer, the starting year for calculation of quantities
#' @param res Positive integer, the temporal resolution (chunks - normally years) over which to calculate quantities
#' @param tsd Character vector of names of types of data: Cat = catch, Ind = relative abundance index, ML = mean length in catches
#' @param stat Character vector of types of quantity to be calculated: slp = slope(log(x)), AAV = average annual variability, mu = mean(log(x))
#' @return A 3D array of results (type of data/stat (e.g. mean catches),time period (chunk), simulation)
#' @author T. Carruthers
#' @references Carruthers and Hordyk 2018
#' @export
getinds<-function(PPD,styr,res=6, tsd= c("Cat","Cat","Cat","Ind","ML"),stat=c("slp","AAV","mu","slp", "slp")){

  nsim<-dim(PPD@Cat)[1]
  proyears<-dim(PPD@Cat)[2]-styr+1

  if(res>proyears)message(paste0("The temporal resolution for posterior predictive data calculation (",res,") is higher than the number of projected years (",proyears,"). Only one time step of indicators are calculated for ",proyears, " projected years."))
  np<-floor(proyears/res)

  ntsd<-length(tsd)
  inds<-array(NA,c(ntsd,np,nsim))

  ind_range <- list()
  for(pp in 1:np) {
    ind<-styr+((pp-1)*res)+1:res
    ind_range[[pp]] <- range(ind)
    for(i in 1:ntsd) inds[i,pp,]<-sapply(1:nsim,get(stat[i]),mat=slot(PPD,tsd[i]),ind=ind)
  }
  #for(i in 1:ntsd){
  #  for(pp in 1:np){
  #    ind<-styr+((pp-1)*res)+1:res
  #    ind_range[[pp]] <- range(ind)
  #    inds[i,pp,]<-sapply(1:nsim,get(stat[i]),mat=slot(PPD,tsd[i]),ind=ind)
  #  }
  #}

  dimnames(inds) <- list(paste0(tsd, "_", stat), vapply(ind_range, paste0, character(1), collapse = "-"), 1:nsim)

  inds

}

#' Produce a cross-correlation plot of the derived data arising from getinds(MSE_object)
#'
#' @param indPPD A 3D array of results arising from running getind on an MSE of the Null operating model (type of data/stat (e.g. mean catches),time period (chunk), simulation)
#' @param indData A 3D array of results arising from running getind on an MSE of the Alternative operating model (type of data/stat (e.g. mean catches),time period (chunk), simulation)
#' @param pp Positive integer, the number of time chunks (blocks of years normally, second dimension of indPPD and indData) to produce the plot for.
#' @param dnam A character vector of names of the data for plotting purposes (as long as dimension 1 of indPPD and indData).
#' @param res The size of the temporal blocking that greated indPPD and indData - this is just used for labelling purposes
#' @return A cross-correlation plot (ndata-1) x (ndata-1)
#' @author T. Carruthers
#' @references Carruthers and Hordyk 2018
#' @export
plot_crosscorr<-function(indPPD,indData,pp=1,dnam=c("CS","CV","CM","IS","MLS"),res=1){

  if(pp>1)namst<-paste(rep(dnam,pp),rep((1:pp)*res,each=length(dnam)))
  if(pp==1)namst=dnam
  cols<-c("#ff000050","#0000ff50")
  ntsd<-dim(indPPD)[1]
  ni<-pp*ntsd
  ind2PPD<-matrix(indPPD[,1:pp,],nrow=ni)
  ind2Data<-matrix(indData[,1:pp,],nrow=ni)
  old_par <- par(no.readonly = TRUE)
  par(mfrow=c(ni-1,ni-1),mai=rep(0,4),omi=c(0.5,0.75,0,0.05))
  on.exit(par(old_par))

  for(i in 2:ni){

    for(j in 1:(ni-1)){

      if(j==i|j>i){

        plot(1,1,col='white',axes=FALSE)

      }else{

        xlim<-quantile(c(ind2PPD[j,],ind2Data[j,]),c(0.02,0.98))
        ylim<-quantile(c(ind2PPD[i,],ind2Data[i,]),c(0.02,0.98))
        plot(ind2PPD[j,],ind2PPD[i,],pch=19,xlim=xlim,ylim=ylim,cex=0.8,col=cols[1],axes=FALSE)
        points(ind2Data[j,],ind2Data[i,],pch=19,cex=0.8,col=cols[2])

      }
      if(i==2&j==(ni-1)){
        legend('center',legend=c("Null OM", "Alternative OM"),text.col=c("blue","red"),bty='n')
      }

      if(j==1)mtext(namst[i],2,line=2,cex=0.6,las=2)
      if(i==ni)mtext(namst[j],1,line=1,cex=0.6,las=2)
      box()

    }

  }

}

#' Calculates mahalanobis distance and rejection of the Null operating model
#'
#' Calculates mahalanobis distance and rejection of the Null operating model, used by wrapping function \link{PRBcalc}.
#'
#' @param indPPD A 3D array of results arising from running getind on an MSE of the Null operating model (type of data/stat (e.g. mean catches),time period (chunk), simulation)
#' @param indData A 3D array of results arising from running getind on an MSE of the Alternative operating model (type of data/stat (e.g. mean catches),time period (chunk), simulation)
#' @param alpha Positive fraction: rate of type I error, alpha
#' @param removedat Logical, should data not contributing to the mahalanobis distance be removed?
#' @param removethresh Positive fraction: the cumulative percentage of removed data (removedat=TRUE) that contribute to the mahalanobis distance
#' @return A list object.
#'  Position 1 is an array of the mahalanobis distances. Dimension 1 is length 2 for the Null OM (indPPD) and the alternative OM (indData).
#'  Dimension 2 is the time block (same length as indPPD dim 2). Dimension 3 is the simulation number (same length at indPPD dim 3.),
#'  Position 2 is a matrix (2 rows, ntimeblock columns) which is (row 1) alpha: the rate of false positives, and row 2 the power (1-beta) the rate of true positives
#' @author T. Carruthers
#' @references Carruthers and Hordyk 2018
#' @importFrom corpcor pseudoinverse
#' @export
Probs<-function(indPPD,indData,alpha=0.05,removedat=FALSE,removethresh=0.05){

  ntsd<-dim(indPPD)[1]
  np<-dim(indPPD)[2]
  nsim<-dim(indPPD)[3]
  PRB<-array(NA,c(2,np))  # False Positive, True Positive
  mah<-array(NA,c(2,np,nsim))

  for(pp in 1:np){

    keep<-array(TRUE,c(ntsd,pp))
    #for(i in 1:ntsd){
    # for(j in 1:pp){
     #  if(dip(indPPD[i,j,])>0.065)keep[i,j]=FALSE
      #}
    #}
    ni<-sum(keep)
    keepind<-as.matrix(expand.grid(1:ntsd,1:pp,1:nsim))[rep(as.vector(keep),nsim),]
    ind3PPD<-t(matrix(indPPD[keepind],nrow=ni))
    ind3Data<-t(matrix(indData[keepind],nrow=ni))

    covr <- cov(ind3PPD)
    mu<-apply(ind3PPD,2,median)

    if(!removedat)keep3<-NA
    if(removedat){
      keep2<-rep(TRUE,ncol(ind3PPD))
      cont<-mahalanobis_contribution(ind3Data,mu,covr)
      mag<-apply(cont,2,mean)
      ind<-order(mag)
      cum<-cumsum(mag[ind])
      keep2[ind[cum<(removethresh*100)]]<-FALSE
      ind3PPD<-ind3PPD[,keep2]
      ind3Data<-ind3Data[,keep2]
      if(pp==1)keep3<-rbind(mag,keep2)
    }

    covr <- cov(ind3PPD)
    mu<-apply(ind3PPD,2,median)

    mahN <- mahalanobis_robust(ind3PPD, center = mu, cov = covr)
    mahA <- mahalanobis_robust(ind3Data, center = mu, cov = covr)

    mahN<-extreme.outlier(mahN)
    mahA<-extreme.outlier(mahA)

    mah[1,pp,]<-mahN
    mah[2,pp,]<-mahA

    Thres<-quantile(mahN,1-alpha,na.rm=T)

    PRB[1,pp]<-mean(mah[1,pp,]>Thres,na.rm=T)   # False positive
    PRB[2,pp]<-mean(mah[2,pp,]>Thres,na.rm=T)   # True positive

  }

  return(list(mahalanobis=mah,PRB=PRB))

}

#' Calculate mahalanobis distance (null and alternative MSEs) and statistical power for all MPs in an MSE
#'
#' @param MSE_null An object of class MSE representing the null hypothesis
#' @param MSE_alt An object of class MSE representing the alternative hypothesis
#' @param tsd Character string of data types: Cat = catch, Ind = relative abundance index, ML = mean length in catches
#' @param stat Character string defining the quantity to be calculated for each data type, slp = slope(log(x)), AAV = average annual variability, mu = mean(log(x))
#' @param dnam Character string of names for the quantities calculated
#' @param res Integer, the resolution (time blocking) for the calculation of PPD
#' @param alpha Probability of incorrectly rejecting the null operating model when it is valid
#' @param plotCC Logical, should the PPD cross correlations be plotted?
#' @param removedat Logical, should data not contributing to the mahalanobis distance be removed?
#' @param removethresh Positive fraction: the cumulative percentage of removed data (removedat=TRUE) that contribute to the mahalanobis distance
#' @importFrom MASS cov.mcd
#' @importFrom corpcor pseudoinverse
#' @return A list object with two hierarchies of indexing, first by MP, second has two positions as described in \link{Probs}: (1) mahalanobis distance, (2) a matrix of type 1 error
#' (first row) and statistical power (second row), by time block.
#' @author T. Carruthers
#' @references Carruthers, T.R, and Hordyk, A.R. In press. Using management strategy evaluation to establish indicators of changing fisheries.
#' Canadian Journal of Fisheries and Aquatic Science.
#' @export
PRBcalc=function(MSE_null,MSE_alt,
                 tsd= c("Cat","Cat","Cat","Ind","ML"),
                 stat=c("slp","AAV","mu","slp", "slp"),
                 dnam=c("C_S","C_V","C_M","I_S","ML_S"),
                 res=6,alpha=0.05,plotCC=FALSE,removedat=FALSE,removethresh=0.025){

  styr<-MSE_null@nyears

  outlist<-new('list')

  for(mm in 1:MSE_null@nMPs){

    if(!is.null(MSE_null@Misc$Data[[mm]])) {
      PPD <- MSE_null@Misc$Data[[mm]]
    } else {
      PPD <- MSE_null@Misc[[mm]]
    }

    if(!is.null(MSE_alt@Misc$Data[[mm]])) {
      Data <- MSE_alt@Misc$Data[[mm]]
    } else {
      Data<-MSE_alt@Misc[[mm]]
    }

    indPPD<-getinds(PPD,styr=styr,res=res,tsd=tsd,stat=stat)
    indData<-getinds(Data,styr=styr,res=res,tsd=tsd,stat=stat)


    if(plotCC) plot_crosscorr(indPPD,indData,pp=2,res=res)

    out<-Probs(indPPD,indData,alpha=alpha,removedat=removedat,removethresh=removethresh)

    dimnames(out[[1]]) <- list(c("Null_OM", "Alt_OM"), dimnames(indPPD)[[2]], dimnames(indPPD)[[3]])
    rownames(out[[2]]) <- c("Type-1", "Power")
    colnames(out[[2]]) <- dimnames(indPPD)[[2]]

    outlist[[mm]]<-out

    message(paste(mm, "of",MSE_null@nMPs,"MPs processed"))

  }

  names(outlist) <- MSE_null@MPs
  return(outlist)

}

mahalanobis_robust<-function (x, center, cov, inverted = FALSE) {

  x <- if (is.vector(x))
    matrix(x, ncol = length(x))
  else as.matrix(x)
  if (!identical(center, FALSE))
    x <- sweep(x, 2L, center)

  invcov <- pseudoinverse(cov)
  setNames(rowSums(x %*% invcov * x), rownames(x))

}


mahalanobis_contribution<-function(ind3Data,mu,covr){

  InvSD <- 1/sqrt(covr[cbind(1:nrow(covr),1:nrow(covr))])
  Dmat<-diag(InvSD)
  DSD <- Dmat %*% covr %*% Dmat
  eig <- eigen(DSD)
  InvRootEig<-1/sqrt(eig$val)
  InvDSDhalf <-  eig$vec %*% diag(InvRootEig) %*% t(eig$vec)
  invcov <- pseudoinverse(covr)
  strcont<-array(NA,c(nrow(ind3Data),ncol(ind3Data)))

  for(i in 1:nrow(ind3Data)){
    diff<-mu-ind3Data[i,]
    Qform<-t(diff)%*% invcov %*% diff
    Tsquare <- Qform[1,1]
    W <- InvDSDhalf %*% Dmat %*% diff
    cont <- diag(W %*% t(W))
    strcont[i,]<-cont/sum(cont)*100
  }

  strcont

}


#' Plot statistical power of the indicator with increasing time blocks
#'
#' @param outlist A list object produced by the function \link{PRBcalc}
#' @param res Integer, the resolution (time blocking) for the calculation of PPD
#' @param maxups Integer, the maximum number of update time blocks to plot
#' @param MPs Character vector of MP names
#' @author T. Carruthers
#' @references Carruthers and Hordyk 2018
#' @export
mahplot<-function(outlist,res=6,maxups=5,MPs){

  nMP<-length(outlist)
  ncol<-floor(nMP^(1/2))
  nrow<-ceiling(nMP/ncol)
  par(mai=c(0.15,0.25,0.25,0.01),omi=c(0.6,0.6,0.01,0.01))
  layout(matrix(1:(nrow*ncol*2),nrow=nrow*2),heights=rep(c(2,1),nrow))
  pmin<-max(0,min(c(0.95,sapply(outlist,function(x)min(x$PRB[2,]))-0.05)))
  plabs<-matrix(paste0("(",letters[1:(nMP*2)],")"),nrow=2)

  for(mm in 1:nMP){
    xaxis=F
    yaxis=F
    if(mm<(nMP/ncol)|mm==(nMP/ncol))yaxis=T
    if(mm%in%((1:ncol)*(nMP/ncol)))xaxis=T

    mahdensplot(outlist[[mm]],xaxis=F,yaxis=yaxis,res=res,maxups=maxups)
    if(yaxis)mtext("Distance D",2,cex=0.9,line=3.5)
    mtext(MPs[mm],3,font=2,line=0.3,cex=0.95)
    mtext(plabs[1,mm],3,line=0.07,adj=0.01,cex=0.88)

    PRBplot(outlist[[mm]],xaxis=xaxis,yaxis=yaxis,res=res,ylim=c(pmin,1),maxups=maxups)
    if(yaxis)mtext("Power",2,cex=0.9,line=3.5)
    mtext(plabs[2,mm],3,line=0.07,adj=0.01,cex=0.88)

  }
  mtext("Future time period",1,outer=T,line=3)

}


mahdensplot<-function(out,adj=0.9,alpha=0.1,xaxis=FALSE,yaxis=FALSE,res=6,maxups=5){

  heightadj<-0.7
  np<-min(dim(out$PRB)[2],maxups)

  densN<-new('list')
  densA<-new('list')
  ylimy=array(0,c(np,2))
  thresh<-rep(0,np)
  tfac<-4
  xref<-rep(0,np)

  xmaxN<-xmaxA<-rep(0,np)


  for(pp in 1:np){
    thresh[pp]<-quantile(out$mah[1,pp,!is.na(out$mah[1,pp,])],1-alpha)
    xref[pp]<-thresh[pp]*tfac
    ylimy[pp,2]<-xref[pp]#quantile(out$mah[1,pp,],0.95,na.rm=T)*1.2
    densN[[pp]]<-density(out$mah[1,pp,],adj=adj,from=0,to=xref[pp],na.rm=T)
    densA[[pp]]<-density(out$mah[2,pp,],adj=adj,from=0,to=xref[pp],na.rm=T)
    ymaxN<-max(densN[[pp]]$y)
    ymaxA<-max(densA[[pp]]$y)
    xmaxN[pp]=ymaxN
    xmaxA[pp]=ymaxA
  }

  for(pp in 1:np){
    #densN[[pp]]$y<-densN[[pp]]$y/max(densN[[pp]]$y)*heightadj
    #densA[[pp]]$y<-densA[[pp]]$y/max(densA[[pp]]$y)*heightadj
    densN[[pp]]$y<-densN[[pp]]$y/xmaxN[pp]*heightadj
    densA[[pp]]$y<-densA[[pp]]$y/xmaxN[pp]*heightadj
  }

  plot(c(0.5,np+0.1),c(0,tfac),col='white',axes=F)

  if(xaxis){
    labs<-paste("Yrs",(1:np)*res-(res-1),"-",(1:np)*res)
    axis(1,1:np,labs)

  }
  if(yaxis){
    axis(2,c(0,1,10E10),c("0","V",""))
  }

  lines(c(-100,100),rep(1,2),lty=2,lwd=1)
  for(pp in 1:np){

    trs<-0.5+pp-1

    lines(densN[[pp]]$y+trs,densN[[pp]]$x/thresh[pp],col='blue')
    lines(densA[[pp]]$y+trs,densA[[pp]]$x/thresh[pp],col='red')

    polyN<-getsegment(densN[[pp]],thresh[pp],lower=T)
    polyA<-getsegment(densA[[pp]],thresh[pp],lower=F)
    polygon(polyN$x+trs,polyN$y/thresh[pp],border=NA,col="#0000ff50")

    if(length(polyA$x)>10)polygon(polyA$x+trs,polyA$y/thresh[pp],border=NA,col="#ff000050")

  }
}


PRBplot<-function(out,res,xaxis=FALSE,yaxis=FALSE,ylim=c(0,1),maxups=5){
  np<-min(dim(out$PRB)[2],maxups)
  plot(c(0.5,np+0.5),range(out$PRB),col='white',axes=F,xlab="",ylab="",ylim=ylim)
  #lines((1:np)-0.2,1-out$PRB[1,],col="#0000ff50",lwd=2)
  abline(h=0.8,lty=2)
  lines((1:np)-0.2,out$PRB[2,1:np],col="#ff000050",lwd=2)
  axis(1,at=c(-1000,10000))
  if(xaxis){
    labs<-paste0("Yrs ",(1:np)*res-(res-1),"-",(1:np)*res)
    at = axTicks(1)
    mtext(side = 1, text = labs, at = at, line = 1,cex=0.9)

    #axis(1,1:np,labs)
  }
  axis(2)
  #if(yaxis)axis(2)

}

getsegment<-function(densobj,thresh,lower=T,inv=T){

  if(lower){
    cond<-densobj$x<thresh
  }else{
    cond<-densobj$x>thresh
  }

  xs<-c(0,densobj$y[cond],0)
  ys<-densobj$x[cond]
  ys<-c(ys[1],ys,ys[length(ys)])

  if(inv)return(list(x=xs,y=ys))
  if(!inv)return(list(x=ys,y=xs))

}

extreme.outlier<-function(x){
  sd<-mean(abs(x-median(x)))
  cond<-(x<(median(x)-3*sd(x)))|(x>(median(x)+3*sd(x)))
  x[cond]<-NA
  x
}

# ====================================================================
