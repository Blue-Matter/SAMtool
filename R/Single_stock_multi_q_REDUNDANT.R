
#' optimize for catchability (q)
#'
#' Function optimizes catchability (q, where F=qE) required to get to user-specified stock
#' depletion
#'
#' @param x Integer, the simulation number
#' @param D A numeric vector nsim long of sampled depletion
#' @param SSB0 A numeric vector nsim long of total unfished spawning biomass
#' @param nareas The number of spatial areas
#' @param maxage The maximum age
#' @param N Array of the numbers-at-age in population. Dimensions are nsim, maxage, nyears, nareas.
#' Only values from the first year (i.e `N[,,1,]`) are used, which is the current N-at-age.
#' @param pyears The number of years to project forward. Equal to 'nyears' for optimizing for q.
#' @param M_ageArray An array (dimensions nsim, maxage, nyears+proyears) with the natural mortality-at-age and year
#' @param Mat_age An array (dimensions nsim, maxage, proyears+nyears) with the proportion mature for each age-class
#' @param Asize A matrix (dimensions nsim, nareas) with size of each area
#' @param Wt_age An array (dimensions nsim, maxage, nyears+proyears) with the weight-at-age and year
#' @param FleetP A list of Fleets containing entries for V, retA, Find and Spat_targ
#' @param mov An array (dimensions nsim, nareas, nareas) with the movement matrix
#' @param SRrel A numeric vector nsim long specifying the recruitment curve to use
#' @param hs A numeric vector nsim long with the steepness values for each simulation
#' @param R0a A matrix (dimensions nsim, nareas) with the unfished recruitment by area
#' @param SSBpR A matrix (dimensions nsim, nareas) with the unfished spawning-per-recruit by area
#' @param aR A numeric vector nareas long with the Ricker SRR a values
#' @param bR A numeric vector nareas long with the Ricker SRR b values
#' @param bounds A numeric vector of length 2 with bounds for the optimizer
#' @param maxF A numeric value specifying the maximum fishing mortality for any single age class
#' @param MPA A matrix of spatial closures by year
#' @param useCPP logical - use the CPP code? For testing purposes only
#' @author T.Carruthers and A. Hordyk
#' @keywords internal
#' @export
getq_multi <- function(x, D, SSB0, nareas, maxage, Np, pyears, M_ageArray, Mat_age, Asize, Wt_age,
                       FleetP, Perr, mov, SRrel, hs, R0a, SSBpR, aR, bR,
                       bounds = c(1e-05, 15),maxF, MPAc, CFp, useCPP=TRUE) {

  # Prep multi-fleet arrays
  nf<-length(FleetP)
  ns<-dim(FleetP[[1]]$V)[1]
  ay<-dim(FleetP[[1]]$V)[3]
  Vuln=retA=array(NA,c(nf,maxage,ay))
  ny<-pyears

  Effind<-t(matrix(unlist(lapply(FleetP,function(dat,x)dat['Find'][[1]][x,],x=x)),ncol=nf))
  Spat_targc<-unlist(lapply(FleetP,function(dat,x)dat['Spat_targ'][[1]][x],x=x))

  for(ff in 1:nf){ # kind of lazy but not THAT slow

    Vuln[ff,,]<-FleetP[[ff]]$V[x,,]
    retA[ff,,]<-FleetP[[ff]]$retA[x,,]

  }

  Fdist<-CFp[x,]/Effind[,ny] # Catch divided by effort (q proxy)
  Fdist<-Fdist/sum(Fdist)    # q ratio proxy (real space)
  par<-c(-5,log(Fdist[2:nf]/(1-Fdist[2:nf]))) # low initial F followed by logit guess at fraction based on Fdist according to catch fraction in recent year

  opt <- optim(par,optQ_multi, method="L-BFGS-B",lower=c(log(bounds[1]),rep(-5,nf-1)),upper=c(log(bounds[2]),rep(5,nf-1)),depc=D[x], SSB0c=SSB0[x], nareas=nareas, maxage=maxage, Ncurr=Np[x,,1,],
               pyears=pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
               Vuln=Vuln, Retc=retA, Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x],
               Effind=Effind,  Spat_targc=Spat_targc, hc=hs[x], R0c=R0a[x,],
               SSBpRc=SSBpR[x,], aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPAc=MPAc, CFc=CFp[x,],useCPP=useCPP,mode='opt')

  out<-optQ_multi(opt$par,depc=D[x], SSB0c=SSB0[x], nareas=nareas, maxage=maxage, Ncurr=Np[x,,1,],
                  pyears=pyears, M_age=M_ageArray[x,,], MatAge=Mat_age[x,,], Asize_c=Asize[x,], WtAge=Wt_age[x,,],
                  Vuln=Vuln, Retc=retA, Prec=Perr[x,], movc=mov[x,,,], SRrelc=SRrel[x],
                  Effind=Effind,  Spat_targc=Spat_targc, hc=hs[x], R0c=R0a[x,], SSBpRc=SSBpR[x,],
                  aRc=aR[x,], bRc=bR[x,], maxF=maxF, MPAc=MPAc, CFc=CFp[x,],useCPP=useCPP,mode='out')

  return(out)

}


#' Optimize q for a single simulation
#'
#' @param logQ log q
#' @param depc Depletion value
#' @param SSB0c Unfished spawning biomass
#' @param nareas Number of areas
#' @param maxage Maximum age
#' @param Ncurr Current N-at-age
#' @param pyears Number of years to project population dynamics
#' @param M_age M-at-age
#' @param Asize_c Numeric vector (length nareas) with size of each area
#' @param MatAge Maturity-at-age
#' @param WtAge Weight-at-age
#' @param Vuln Vulnerability-at-age
#' @param Retc Retention-at-age
#' @param Prec Recruitment error by year
#' @param movc movement matrix
#' @param SRrelc SR parameter
#' @param Effind Historical fishing effort
#' @param Spat_targc Spatial targetting
#' @param hc Steepness
#' @param R0c Unfished recruitment by area
#' @param SSBpRc Unfished spawning biomass per recruit by area
#' @param aRc Ricker aR
#' @param bRc Ricker bR
#' @param maxF maximum F
#' @param MPA A matrix of spatial closures by year
#' @param useCPP Logical. Use the CPP code?
#' @param model Character 'opt' is optimization and returns an objective funciton to be minimized, anything else returns fitted values
#' @author A. Hordyk
#' @keywords internal
#' @export
optQ_multi <- function(par, depc, SSB0c, nareas, maxage, Ncurr, pyears, M_age,
                       MatAge, Asize_c, WtAge, Vuln, Retc, Prec, movc, SRrelc, Effind, Spat_targc, hc,
                       R0c, SSBpRc, aRc, bRc, maxF, MPAc, CFc, useCPP, mode = 'opt') {


  nf<-dim(Effind)[1]
  ny<-dim(Effind)[2]
  ay<-dim(Vuln)[3]

  qtot<-exp(par[1])
  qlogit<-c(0,par[2:nf])
  qfrac<-exp(qlogit)/(sum(exp(qlogit)))
  Effdist<-qfrac*Effind
  Efftot<-apply(Effdist,2,sum)

  MPAind<-TEG(c(nf,ny,nareas))
  MPAtemp<-array(1/nf,dim(MPAc))
  MPAtemp[MPAind]=(MPAc[MPAind]*Effdist[MPAind[,1:2]])/Efftot[MPAind[,2]] # weighted by effort and fleet exposure
  MPAf<-apply(MPAtemp,2:3,sum)

  Vulnf<-Retf<-array(NA,c(maxage,ay))
  Vind<-TEG(c(nf,maxage,ny))
  Vtemp<-Vuln[,,1:ny]
  Vtemp[Vind]<-(Vuln[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
  Vulnf[,1:ny]<-apply(Vtemp,2:3,sum)
  Vulnf[,(ny+1):ay]<-Vulnf[,ny]  # Future vulnerability is the same
  Vulnf<-nlz(Vulnf,2,"max")       # normalize to max 1

  Rtemp<-Retc[,,1:ny]
  Rtemp[Vind]<-(Retc[Vind]*Effdist[Vind[,c(1,3)]])/Efftot[Vind[,3]]
  Retf[,1:ny]<-apply(Rtemp,2:3,sum)
  Retf[,(ny+1):ay]<-Retf[,ny]   # Future retention is the same

  Spat_targf<-sum(apply(Effdist,1,sum)*Spat_targc)/sum(Effdist) # Approximation according to historical F by fleet

  simpop <- DLMtool::popdynCPP(nareas, maxage, Ncurr, pyears, M_age, Asize_c,
                      MatAge, WtAge, Vuln=Vulnf, Retc=Retf, Prec, movc, SRrelc, Effind=Efftot, Spat_targc=Spat_targf, hc,
                      R0c=R0c, SSBpRc=SSBpRc, aRc=aRc, bRc=bRc, Qc=qtot, Fapic=0,
                      maxF=maxF, MPA=MPAf, control=1,  SSB0c=SSB0c)



  # 1:Narray 2:Barray 3:SSNarra 4:SBarray 5:VBarray 6:FMarray 7:FMretarray 8:Zarray
  ssb <- sum(simpop[[4]][,pyears,])
  BB<- apply(simpop[[2]][,pyears,],1,sum)  # age area
  Cf<-array(rep(BB,each=nf)*Vuln[,,ny]*Effdist[,ny],c(nf,maxage))
  Cpred<-apply(Cf,1,sum)/sum(Cf)

  depOBJ<-50*(log(depc) - log(ssb/SSB0c))^2
  cOBJ<-sum(log(Cpred[2:nf]/CFc[2:nf])^2)

  if(mode=='opt'){
    return(depOBJ+cOBJ)
  }else{
    return(list(qtot=qtot,qfrac=qfrac,CFc=CFc,Cpred=Cpred,depc=depc,deppred=ssb/SSB0c,Vulnf=Vulnf,Retf=Retf,MPAf=MPAf))
  }

}


