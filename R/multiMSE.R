# =================================================================================================
# ====== multiMSE  - an R function for doing MSE evaluation for a multi-stock multi-fleet OM ======
# =================================================================================================

#' Run a multi-fleet multi-stock Management Strategy Evaluation
#'
#' A function that runs a Management Strategy Evaluation (closed-loop
#' simulation) for a specified operating model
#'
#'
#' @param MOM A multi-fleet multi-stock operating model (class 'MOM')
#' @param MPs A matrix of methods (nstock x nfleet) (character string) of class MP
#' @param CheckMPs Logical to indicate if \link{Can} function should be used to check
#' if MPs can be run.
#' @param timelimit Maximum time taken for a method to carry out 10 reps
#' (methods are ignored that take longer)
#' @param Hist Should model stop after historical simulations? Returns a list
#' containing all historical data
#' @param ntrials Maximum of times depletion and recruitment deviations are
#' resampled to optimize for depletion. After this the model stops if more than
#' percent of simulations are not close to the required depletion
#' @param fracD Maximum allowed proportion of simulations where depletion is not
#' close to sampled depletion from OM before model stops with error
#' @param CalcBlow Should low biomass be calculated where this is the spawning
#' biomass at which it takes HZN mean generation times of zero fishing to reach
#' Bfrac fraction of SSBMSY
#' @param HZN The number of mean generation times required to reach Bfrac SSBMSY
#' in the Blow calculation
#' @param Bfrac The target fraction of SSBMSY for calculating Blow
#' @param AnnualMSY Logical. Should MSY statistics be calculated for each projection year?
#' May differ from MSY statistics from last historical year if there are changes in productivity
#' @param silent Should messages be printed out to the console?
#' @param PPD Logical. Should posterior predicted data be included in the MSE object Misc slot?
#' @param parallel Logical. Should the MSE be run using parallel processing?
#' @param save_name Character. Optional name to save parallel MSE list
#' @param checks Logical. Run tests?
#' @param control control options for testing and debugging
#' @return A hierarchical list (by stock then fleet) of objects of class \linkS4class{MSE}
#' @author T. Carruthers and A. Hordyk
#' @export
multiMSE <- function(MOM, MPs = list(c("AvC","DCAC"),c("FMSYref","curE")),
                   CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=FALSE,
                   HZN=2, Bfrac=0.5, AnnualMSY=TRUE, silent=FALSE, PPD=FALSE, parallel=FALSE,
                   save_name=NULL, checks=FALSE, control=NULL) {

  # Known issues:
  # 1) maxage is the same across species so computationally inefficient
  # 2) initdist is no yet working as a cpars value (indeed I'm not sure where this comes into the OM specification at all!)
  # 3) NO custom MP check whether already in DLMtool
  # 4) No parallel processing option at present

  if (class(MOM)!='MOM') stop("OM object is not of class '<OM'", call. = FALSE)
  #if (Hist & parallel) {
  #  message("Sorry! Historical simulations currently can't use parallel.")
   # parallel <- FALSE
  #}

  if (parallel) {
    message("Sorry! multiMSE currently can't use parallel.")
    parallel <- FALSE
  }

  # Set DLMenv to be empty. Currently updated by Assess models in MSEtool
  rm(list = ls(DLMenv), envir = DLMenv)

  # Check MPs
  MPvec<-unlist(MPs)
  if (!all(is.na(MPvec))) {
    for (mm in MPvec) {
      chkMP <- try(get(mm), silent=TRUE)
      if (!(class(chkMP) %in% c('MP','MMP'))) stop(mm, " is not a valid MP", call.=FALSE)
    }
  }


  MSE1 <- multiMSE_int(MOM, MPs, CheckMPs, timelimit, Hist, ntrials, fracD, CalcBlow,
                       HZN, Bfrac, AnnualMSY, silent, PPD, checks=checks, control=control)


  return(MSE1)

}

multiMSE_int <- function(MOM, MPs=list(c("AvC","DCAC"),c("FMSYref","curE")),
                       CheckMPs = FALSE, timelimit = 1, Hist=FALSE, ntrials=50, fracD=0.05, CalcBlow=FALSE,
                       HZN=2, Bfrac=0.5, AnnualMSY=TRUE, silent=FALSE, PPD=FALSE, checks=FALSE,
                       control=NULL, parallel=FALSE) {

  # Priority list
  # Multi-fleet popdynCPP needed to properly account for fleet-specific MPAs, spatial targetting etc.

  # Known issues
  # 1) SSN[SPAYR] <- Nfrac[SA] * StockPars[[p]]$R0[S] * StockPars[[p]]$initdist[SAR]  this is Pinitdist[SR] in runMSE but I have no idea where Pinitdist comes from
  # 2) currently assume unfished vulnerability is equally weighted among fleets (first V calculation)
  # 3) MICE mode: MSY calcs are fully dynamic and by year - we do not solve the long-term equilibrium MSY for MICE models
  # 4) No check for MOM correct object formatting yet
  # 5) MICE mode: SSB0 is constant (does not change with M for example) and based on long-term ecosystem average as specified in StockPars[[p]]$SSB0
  # 6) For MSY calculations, vulnerability is calculated by fishing mortality rate summed over both areas (not weighted by total catches for example)
  # 7) Blow calculations are currently not coded!
  # 8) RefY reference yield is currently MSY!
  # 9) No control$Cbias_yr functionality!
  # 10) MGT calculation not currently coded!
  # 11) No annual MSY calculation currently (need updating inside MP loop because it now responds to MICE modelled relationships)
  # 12) Currently there is no MMP checking
  # 13) SelectChanged, the switch that only does MSY calcs if selectivity has changed, is currently set permanently to TRUE (around line 1392)
  # 14) Historical vulnerability is calculated across fleets according to todays catch fraction
  # 15) Pstar quantiles of TACs not available nreps fixed to 1

  # Needs checking
  # 1) The vulnerability in the fleets does not max to 1
  # 2) When calculating aggregate retention in the multi_q_estimation function, how should max retention at age be calculated (what is the max value)?
  # 3) Single simulation run nsim=1
  # 4) get rid of any DLMtool::: triple colons

  if (class(MOM) != "MOM") stop("You must specify a valid operating model of class MOM (multi operating model)")
  Misc<-new('list') # Blank miscellaneous slot created
  if("seed"%in%slotNames(MOM)) set.seed(MOM@seed) # set seed for reproducibility
  if(is.na(MPs[1]))stop("You must specify valid management procedures - argument MPs")

  #OM <- updateMSE(OM)
  tiny <- 1e-15  # define tiny variable
  nsim<-MOM@nsim
  nyears <- MOM@Fleets[[1]][[1]]@nyears  # number of historical years
  proyears<-MOM@proyears
  allyears<-nyears+proyears
  Stocks<-MOM@Stocks
  Fleets<-MOM@Fleets
  Obs<-MOM@Obs
  Imps<-MOM@Imps
  Rel<-MOM@Rel
  SexPars<-MOM@SexPars
  Complexes<-MOM@Complexes
  CatchFrac<-MOM@CatchFrac

  np<-length(Stocks)
  nf<-length(Fleets[[1]])

  if(np==1&nf==1){
    message("You have specified only a single stock and fleet. You should really be using the function runMSE()")
  }else if(np>1 & length(MOM@Rel)==0 & length(MOM@SexPars)==0){
    message("You have specified more than one stock but no MICE relationships (slot MOM@Rel) or sex-specific relationships (slot MOM@SexPars) among these. As they are independent, consider doing MSE for one stock at a time for computational efficiency.")
  }

  maxF<-MOM@maxF
  Snames<-SIL(Stocks,"Name")
  Fnames<-matrix(make.unique(SIL(Fleets,"Name")),nrow=nf)
  cpars<-MOM@cpars

  #MOM <- ChkObj(MOM) # Check that all required slots in OM object contain values
  if (proyears < 2) stop('OM@proyears must be > 1', call.=FALSE)
  ### Sampling OM parameters ###
  if(!silent) message("Loading operating model")

  # Allocation
  if(length(MOM@Allocation)==0){
    MOM@Allocation <-CatchFrac
    message("Slot @Allocation of MOM object not specified. Setting slot @Allocation equal to slot @CatchFrac - current catch fractions")
  }

  if(length(MOM@Efactor)==0){
    MOM@Efactor <-list()
    for(p in 1:np)MOM@Efactor[[p]]<-array(1,c(nsim,nf))
    message("Slot @Efactor of MOM object not specified. Setting slot @Efactor to current effort for all fleets")
  }

  # --- Sample custom parameters ----

  SampCpars<-list() # empty list

  # custom parameters exist - sample and write to list
  for(p in 1:np){
    SampCpars[[p]]<-list()
    for(f in 1:nf){
      if(length(cpars[[p]][[f]])>0){
        message(paste(Stocks[[p]]@Name," - ",Fleets[[p]][[f]]@Name))
        ncparsim<-cparscheck(cpars[[p]][[f]])   # check each list object has the same length and if not stop and error report
        SampCpars[[p]][[f]] <- SampleCpars(cpars[[p]][[f]], nsim, msg=!silent)
      }else{
        SampCpars[[p]][[f]] <-list()
      }
    }
  }

  # All stocks and sampled parameters must have compatible array sizes (maxage)
  maxage_s<-unique(SIL(MOM@Stocks,"maxage"))
  if(length(maxage_s)>1)message(paste("Stocks of varying maximum ages have been specified, all simulations will run to",max(maxage_s),"ages"))
  maxage<-max(maxage_s)
  for(p in 1:np)MOM@Stocks[[p]]@maxage<-maxage

  # --- Sample Stock Parameters ----
  StockPars<-FleetPars<-ObsPars<-ImpPars<-new('list')
  for(p in 1:np){
    StockPars[[p]] <- SampleStockPars(MOM@Stocks[[p]], nsim, nyears, proyears,  SampCpars[[p]][[1]], msg=!silent)
  }

  # --- Sample Fleet Parameters ----
  for(p in 1:np){
    FleetPars[[p]]<-ObsPars[[p]]<-ImpPars[[p]]<-list()
    for(f in 1:nf){
      FleetPars[[p]][[f]] <- SampleFleetPars(MOM@Fleets[[p]][[f]], Stock=StockPars[[p]], nsim, nyears, proyears, cpars=SampCpars[[p]][[f]])
      #FleetPars[[p]][[f]]$Find<-FleetPars[[p]][[f]]$Find/(apply(FleetPars[[p]][[f]]$Find,1,mean))
    }
  }
  for(p in 1:np)for(f in 1:nf)ObsPars[[p]][[f]] <- SampleObsPars(MOM@Obs[[p]][[f]], nsim, cpars=SampCpars[[p]][[f]])
  for(p in 1:np)for(f in 1:nf) ImpPars[[p]][[f]] <- SampleImpPars(MOM@Imps[[p]][[f]],nsim, cpars=SampCpars[[p]][[f]])


  # Assign Stock pars to function environment
  #for (X in 1:length(StockPars)) assign(names(StockPars)[X], StockPars[[X]])

  # Assign Fleet pars to function environment
  #for (X in 1:length(FleetPars)) assign(names(FleetPars)[X], FleetPars[[X]])

  # Assign Obs pars to function environment
  # for (X in 1:length(ObsPars)) assign(names(ObsPars)[X], ObsPars[[X]])

  # Assign Imp pars to function environment
  # for (X in 1:length(ImpPars)) assign(names(ImpPars)[X], ImpPars[[X]])

  # Check for compatibility of dimensions among stocks

  nareas_s<-NIL(StockPars,"nareas",lev1=T)
  if(length(unique(nareas_s))!=1)stop("Stocks must have the same specified number of areas - check cpars$mov for each stock object")
   nareas<-nareas_s[1]

  ### End of sampling OM parameters ###
  # nsim, ns, nf, maxage, nyears, nareas
  N <- Biomass <- Z<- VBiomass<- SSN <- SSB <- array(NA, dim = c(nsim, np, maxage, nyears, nareas))
  VF<-FretA<-array(NA, dim = c(nsim, np, nf, maxage, allyears))
  VBF<-FM <- FMret <- array(NA, dim = c(nsim, np, nf, maxage, nyears, nareas))  # stock numbers array
  SPR <- array(NA, dim = c(nsim, np, maxage, nyears)) # store the Spawning Potential Ratio
  MPA<-array(1,c(np,nf, nyears+proyears,nareas))
  Agearray <- array(rep(1:maxage, each = nsim), dim = c(nsim, maxage))  # Age array

  # Hermaphroditism (this is the fraction to be kept (after sex change)) E.g. protygynous (Female to male) is H_1_2 where 1 is female 2 is male
  HermFrac<-expandHerm(SexPars$Herm,maxage=StockPars[[1]]$maxage,np=np,nsim=nsim) # [sim, stock, maxage] Defaults to all 1s if length(SexPars$Herm)==0

  for(p in 1:np){

    surv <- matrix(1, nsim, maxage)
    surv[, 2:maxage] <- t(exp(-apply(StockPars[[p]]$M_ageArray[,,1], 1, cumsum)))[, 1:(maxage-1)]  # Survival array

    Nfrac <- surv * StockPars[[p]]$Mat_age[,,1] * HermFrac[,p,]  # predicted Numbers of mature ages in first year

    SPAYR <- as.matrix(expand.grid(1:nareas, 1, 1:maxage, p, 1:nsim)[5:1])  # Set up some array indexes sim (S) age (A) year (Y) region/area (R)
    SPA<-SPAYR[,1:3]
    SAY <- SPAYR[, c(1,3,4)]
    SAR <- SPAYR[, c(1,3,5)]
    SA <- Sa <- SPAYR[, c(1,3)]
    SR <- SPAYR[, c(1,5)]
    S <- SPAYR[, 1]
    SY <- SPAYR[, c(1, 4)]
    Sa[,2]<-maxage-Sa[,2]+1 # This is the process error index for initial year

    # <<<< Previously an if statement for existence of initdist, see runMSE.R of DLMtool/R>>>>>>>>

    R0a <- matrix(StockPars[[p]]$R0, nrow=nsim, ncol=nareas, byrow=FALSE) * StockPars[[p]]$initdist[,1,] #*HermFrac[,p,1]  # !!!! INITDIST OF AGE 1. Unfished recruitment by area

    SSN[SPAYR] <- Nfrac[SA] * StockPars[[p]]$R0[S] * StockPars[[p]]$initdist[SAR]  # Calculate initial spawning stock numbers
    N[SPAYR] <- StockPars[[p]]$R0[S] * surv[SA] *HermFrac[SPA] * StockPars[[p]]$initdist[SAR]  # Calculate initial stock numbers
    Neq <- N
    Biomass[SPAYR] <- N[SPAYR] * StockPars[[p]]$Wt_age[SAY]  # Calculate initial stock biomass
    SSB[SPAYR] <- SSN[SPAYR] * StockPars[[p]]$Wt_age[SAY]    # Calculate spawning stock biomass

    Vraw<-array(NIL(listy=FleetPars[[p]],namey="V"),c(nsim,maxage,allyears,nf))
    Vind<-as.matrix(expand.grid(1:nsim,p,1:nf,1:maxage,1:allyears))
    VF[Vind]<-Vraw[Vind[,c(1,4,5,3)]]

    if(nf==1){
      V<-VF[,p,1,,] #<-SOL(FleetPars[[p]],"V")
    }else{
      #Weight by catch fraction
      V<-array(0,c(nsim,maxage,allyears))
      for(f in 1:nf){
        V<-V+VF[,p,f,,]*CatchFrac[[p]][,f]
      }
      #V<-nlz(V,c(1,3),"max") # currently assume unfished vulnerability is equally weighted among fleets
      # V includes discards
    }

    VBiomass[SPAYR] <- Biomass[SPAYR] * V[SAY]  # Calculate vunerable biomass

    Fretraw<-array(NIL(listy=FleetPars[[p]],namey="retA"),c(nsim,maxage,allyears,nf))
    FretA[Vind]<-Fretraw[Vind[,c(1,4,5,3)]]

    if (nsim > 1) {
      SSN0 <- apply(SSN[,p , , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
      SSB0 <- apply(SSB[,p , , 1, ], 1, sum)  # Calculate unfished spawning stock biomass
      SSBpR <- matrix(SSB0/StockPars[[p]]$R0, nrow=nsim, ncol=nareas)  # Spawning stock biomass per recruit
      SSB0a <- apply(SSB[,p, , 1, ], c(1, 3), sum)  # Calculate unfished spawning stock numbers
      B0 <- apply(Biomass[,p, , 1, ], 1, sum)
      N0 <- apply(N[,p, , 1, ], 1, sum)
    } else {
      SSN0 <- apply(SSN[,p, , 1, ], 2, sum)  # Calculate unfished spawning stock numbers
      SSB0 <-  sum(SSB[,p, , 1, ])  # Calculate unfished spawning stock biomass
      SSBpR <- SSB0/R0  # Spawning stock biomass per recruit
      SSB0a <- apply(SSB[,p, , 1, ], 2, sum)  # Calculate unfished spawning stock numbers
      B0 <- apply(Biomass[,p, , 1, ], 2, sum)
      N0 <- apply(N[,p, , 1, ], 2, sum)
    }

    bR <- matrix(log(5 * StockPars[[p]]$hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
    aR <- matrix(exp(bR * SSB0a)/SSBpR, nrow=nsim)  # Ricker SR params

    #  --- Non-equilibrium calcs ----
    SSN[SPAYR] <- Nfrac[SA] * StockPars[[p]]$R0[S] * StockPars[[p]]$initdist[SAR]*StockPars[[p]]$Perr_y[Sa]  # Calculate initial spawning stock numbers
    N[SPAYR] <- StockPars[[p]]$R0[S] * surv[SA] * HermFrac[SPA]* StockPars[[p]]$initdist[SAR]*StockPars[[p]]$Perr_y[Sa]  # Calculate initial stock numbers

    Biomass[SPAYR] <- N[SPAYR] * StockPars[[p]]$Wt_age[SAY]  # Calculate initial stock biomass
    SSB[SPAYR] <- SSN[SPAYR] * StockPars[[p]]$Wt_age[SAY]    # Calculate spawning stock biomass
    VBiomass[SPAYR] <- Biomass[SPAYR] * V[SAY]  # Calculate vunerable biomass

    # Assign stock parameters to StockPars object

    StockPars[[p]]$SSBpR <-SSBpR
    StockPars[[p]]$aR <-aR
    StockPars[[p]]$bR<-bR
    StockPars[[p]]$SSB0<-SSB0
    StockPars[[p]]$R0a <-R0a
    StockPars[[p]]$surv<-surv
    StockPars[[p]]$B0<-B0
    StockPars[[p]]$N0<-N0

    for(f in 1:nf)FleetPars[[p]][[f]]$V<-VF[,p,f,,]

    # --- Historical Spatial closures ----

    for(f in 1:nf){

      if (all(!is.na(Fleets[[p]][[f]]@MPA)) && sum(Fleets[[p]][[f]]@MPA) != 0) {

        yrindex <- Fleets[[p]][[f]]@MPA[,1]
        if (max(yrindex)>nyears) stop("Invalid year index for spatial closures: must be <= nyears")
        if (min(yrindex)<1) stop("Invalid year index for spatial closures: must be > 1")
        if (ncol(Fleets[[p]][[f]]@MPA)-1 != nareas) stop("OM@MPA must be nareas + 1")

        for (xx in seq_along(yrindex)) {
          MPA[p,f,yrindex[xx]:nrow(MPA),] <- matrix(Fleets[[p]][[f]]@MPA[xx, 2:ncol(Fleets[[p]][[f]]@MPA)], nrow=length(yrindex[xx]:nrow(MPA)),ncol=nareas, byrow = TRUE)
        }

      }
    }

  } # end of np


  # if SexPars Sex specific exceptions (recalculation of SSB0, aR, bR, SSBpR)
  if(length(SexPars)>0){

    message("You have specified sex-specific dynamics, unfished spawning biomass and specified stock depletion will be mirrored across sex types according to SexPars$SSBfrom")

    SSB0s<-matrix(NIL(StockPars,"SSB0"),nrow=nsim) # sim, p
    sexmatches<-sapply(1:nrow(SexPars$SSBfrom),function(x,mat)paste(mat[x,],collapse="_"), mat=SexPars$SSBfrom)
    parcopy<-match(sexmatches,sexmatches)
    StockPars_t<-StockPars # need to store a temporary object for copying to/from

    for(p in 1:np){

      SSB0<-apply(matrix(rep(SexPars$SSBfrom[p,],each=nsim),nrow=nsim)*SSB0s,1,sum)
      StockPars[[p]]$SSB0<-SSB0

      # copied parameters
      StockPars[[p]]$D<-StockPars_t[[parcopy[p]]]$D
      StockPars[[p]]$hs<-StockPars_t[[parcopy[p]]]$hs
      StockPars[[p]]$R0<-StockPars_t[[parcopy[p]]]$R0
      StockPars[[p]]$R0a<-StockPars_t[[parcopy[p]]]$R0a

      StockPars[[p]]$SSBpR<-array(SSB0/StockPars[[p]]$R0,c(nsim,nareas)) # !!!!!!!!!!! SSBpR hardwired to be the same among areas !!!!

      idist<-StockPars[[p]]$R0a/apply(StockPars[[p]]$R0a,1,sum)
      SSB0a<-SSB0*idist

      StockPars[[p]]$bR <- matrix(log(5 * StockPars[[p]]$hs)/(0.8 * SSB0a), nrow=nsim)  # Ricker SR params
      StockPars[[p]]$aR <- matrix(exp(StockPars[[p]]$bR * SSB0a)/StockPars[[p]]$SSBpR, nrow=nsim)  # Ricker SR params


    }

    if(length(SexPars$Herm)>0){
      message("You have specified sequential hermaphroditism (SexPars$Herm). Unfished stock numbers will be calculated from this vector of fractions at age. Population dynamics will move individuals from one sex to another")
    }

  }

  # --- Optimize catchability (q) to fit depletion ----
  #if ('unfished' %in% names(control) && control$unfished) {
  # if(!silent) message("Simulating unfished historical period")
  #  Hist <- TRUE
  # CalcBlow <- FALSE
  #  qs <- array(0,dim(nsim,nf)) # no fishing
  #} else {

  if(!silent) message("Optimizing for user-specified depletion (takes approximately [(nstocks x nfleets)/(9 x number of cores in cluster)] minutes per simulation)")  # Print a progress update

  bounds <- c(0.0001, 15) # q bounds for optimizer

  if(snowfall::sfIsRunning()){
    out<-snowfall::sfLapply(1:nsim,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF,
                            MPA,CatchFrac, bounds= bounds,tol=1E-6,Rel,SexPars)
  }else{
    out<-lapply(1:nsim,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF,
                MPA,CatchFrac, bounds= bounds,tol=1E-6,Rel,SexPars)
  }

  # system.time({ out<-lapply(1:1,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF, MPA,CatchFrac, bounds= bounds,tol=1E-6,Rel)})
  # system.time({out<-sfLapply(1:nsim,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF, MPA,CatchFrac, bounds= bounds,tol=1E-6,Rel)})

  qs<-t(matrix(NIL(out,"qtot"),nrow=np))
  qfrac<-aperm(array(NIL(out,"qfrac"),c(np,nf,nsim)),c(3,1,2))

  for(p in 1:np){
    for(f in 1:nf){
      FleetPars[[p]][[f]]$qs<-qs[,p]*qfrac[,p,f]
    }
  }

  # --- Check that q optimizer has converged ----
  LimBound <- c(1.1, 0.9)*range(bounds)  # bounds for q (catchability). Flag if bounded optimizer hits the bounds
  probQ <- which(apply(qs > max(LimBound) | qs < min(LimBound),1,sum)>0)
  Nprob <- length(probQ)

  # If q has hit bound, re-sample depletion and try again. Tries 'ntrials' times and then alerts user
  if (length(probQ) > 0) {

    Err <- TRUE
    if(!silent) message(Nprob,' simulations have final biomass that is not close to sampled depletion')
    if(!silent) message('Re-sampling depletion, recruitment error, and fishing effort')

    count <- 0
    MOM2 <- MOM

    while (Err & count < ntrials) {
      # Re-sample Stock Parameters

      Nprob <- length(probQ)
      MOM2@nsim <- Nprob

      SampCpars2 <- vector("list", nf)

      for(p in 1:np){

        for(f in 1:nf){
          if(length(cpars[[p]][[f]])>0){
            #message(paste(Stocks[[p]]@Name," - ",Fleets[[p]][[f]]@Name))
            ncparsim<-cparscheck(cpars[[p]][[f]])   # check each list object has the same length and if not stop and error report
            SampCpars2[[f]] <- SampleCpars(cpars[[p]][[f]], Nprob, msg=!silent)
          }
        }

        #StockPars[[p]] <- SampleStockPars(MOM@Stocks[[p]], nsim, nyears, proyears,  SampCpars[[p]][[f]], msg=!silent)
        ResampStockPars <- SampleStockPars(MOM2@Stocks[[p]], nsim=Nprob,nyears=nyears,proyears=proyears,cpars=SampCpars2[[1]], msg=FALSE)
        #ResampStockPars$CAL_bins <- StockPars$CAL_bins
        #ResampStockPars$CAL_binsmid <- StockPars$CAL_binsmid

        # Re-sample depletion
        StockPars[[p]]$D[probQ] <- ResampStockPars$D

        # Re-sample recruitment deviations
        StockPars[[p]]$procsd[probQ] <- ResampStockPars$procsd
        StockPars[[p]]$AC[probQ] <- ResampStockPars$AC
        StockPars[[p]]$Perr_y[probQ,] <- ResampStockPars$Perr_y
        StockPars[[p]]$hs[probQ] <- ResampStockPars$hs
      } # end of P
      # Re-sample historical fishing effort

      ResampFleetPars<- vector("list", nf)
      for(p in 1:np){
        for(f in 1:nf){
          ResampFleetPars <- SampleFleetPars(MOM2@Fleets[[p]][[f]], Stock=ResampStockPars, nsim=Nprob, nyears=nyears, proyears=proyears, cpars=SampCpars2[[f]])
          FleetPars[[p]][[f]]$Esd[probQ] <- ResampFleetPars$Esd
          FleetPars[[p]][[f]]$Find[probQ, ] <- ResampFleetPars$Find
          FleetPars[[p]][[f]]$dFfinal[probQ] <- ResampFleetPars$dFfinal
        }
      }

      if(snowfall::sfIsRunning()){
        out2<-snowfall::sfLapply(probQ,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF,
                                MPA,CatchFrac, bounds= bounds,tol=1E-6,Rel,SexPars)
      }else{
        out2<-lapply(probQ,getq_multi_MICE,StockPars, FleetPars, np,nf, nareas, maxage, nyears, N, VF, FretA, maxF=MOM@maxF,
                    MPA,CatchFrac, bounds= bounds,tol=1E-6,Rel,SexPars)
      }

      qs2<-t(matrix(NIL(out2,"qtot"),nrow=np))
      qout2<-array(NIL(out2,"qfrac"),c(np,nf,nsim))
      qfrac2<-array(NA,c(Nprob,np,nf))
      qind2<-TEG(dim(qfrac2))
      qfrac2[qind2]<-qout2[qind2[,c(2,3,1)]]
      qfrac[probQ,,]<-qfrac2
      qs[probQ,]<-qs2

      probQ <- which(apply(qs > max(LimBound) | qs < min(LimBound),1,sum)>0)
      count <- count + 1
      if (length(probQ) == 0) Err <- FALSE

    }
    if (Err) { # still a problem

      tooLow <- length(which(qs > max(LimBound)))
      tooHigh <- length(which(qs < min(LimBound)))
      prErr <- length(probQ)/nsim
      if (prErr > fracD & length(probQ) >= 1) {
        if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
        if(!silent) message("More than ", fracD*100, "% of simulations can't get to the specified level of depletion with these Operating Model parameters")
        stop("Change OM@seed and try again for a complete new sample, modify the input parameters, or increase ntrials")
      } else {
        if (length(tooLow) > 0) message(tooLow, " sims can't get down to the lower bound on depletion")
        if (length(tooHigh) > 0) message(tooHigh, " sims can't get to the upper bound on depletion")
        if(!silent) message("More than ", 100-fracD*100, "% simulations can get to the sampled depletion.\nContinuing")
      }
    }

    for(p in 1:np)for(f in 1:nf) FleetPars[[p]][[f]]$qs<-qs[,p]*qfrac[,p,f]

  }  # end of q estimation steps

  if(!silent) message("Calculating historical stock and fishing dynamics")  # Print a progress update

  histYrs <- sapply(1:nsim,HistMICE, StockPars=StockPars,FleetPars=FleetPars,np=np,nf=nf,nareas=nareas,
                    maxage=maxage,nyears=nyears,N=N,VF=VF,FretA=FretA,maxF=MOM@maxF,MPA=MPA,Rel=Rel,SexPars=SexPars,qs=qs,qfrac=qfrac)

  N <- aperm(array(as.numeric(unlist(histYrs[1,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))
  Biomass <- aperm(array(as.numeric(unlist(histYrs[2,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))
  SSN <- aperm(array(as.numeric(unlist(histYrs[3,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))
  SSB <- aperm(array(as.numeric(unlist(histYrs[4,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))
  VBiomass <- aperm(array(as.numeric(unlist(histYrs[5,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))
  FM <- aperm(array(as.numeric(unlist(histYrs[6,], use.names=FALSE)), dim=c(np ,nf,maxage, nyears, nareas, nsim)), c(6,1,2,3,4,5))
  FMret <- aperm(array(as.numeric(unlist(histYrs[7,], use.names=FALSE)), dim=c(np ,nf,maxage, nyears, nareas, nsim)), c(6,1,2,3,4,5))
  VBF <- aperm(array(as.numeric(unlist(histYrs[15,], use.names=FALSE)), dim=c(np ,nf,maxage, nyears, nareas, nsim)), c(6,1,2,3,4,5))
  Z <- aperm(array(as.numeric(unlist(histYrs[16,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))
  FMt<-aperm(array(as.numeric(unlist(histYrs[17,], use.names=FALSE)), dim=c(np ,maxage, nyears, nareas, nsim)), c(5,1,2,3,4))

  # Depletion check
  SSB0_specified <- array(NIL(StockPars,'SSB0'),c(nsim,np))
  D_specified <- array(NIL(StockPars,'D'),c(nsim,np))
  Depletion <- apply(SSB[,,,nyears,,drop=F],1:2,sum)/ SSB0_specified

  if(length(SexPars)>0){ # need to copy over depletion for a sex-specific model

    sexmatches<-sapply(1:nrow(SexPars$SSBfrom),function(x,mat)paste(mat[x,],collapse="_"), mat=SexPars$SSBfrom)
    parcopy<-match(sexmatches,sexmatches)
    StockPars_t<-StockPars # need to store a temporary object for copying to/from
    Depletion[,1:np]<-Depletion[,parcopy]

  }


  # if (nsim == 1) Depletion <- sum(SSB[,p,,nyears,])/SSB0 #^betas
  for(p in 1:np)StockPars[[p]]$Depletion<-Depletion[,p]

  if(checks){
    Cpred<-array(NA,c(nsim,np,nf,maxage,nareas))
    Cind<-as.matrix(expand.grid(1:nsim,1:np,1:nf,1:maxage,nyears,1:nareas))
    Cpred[Cind[,c(1:4,6)]]<-Biomass[Cind[,c(1,2,4,5,6)]]*(1-exp(-FM[Cind]))
    Cpred<-apply(Cpred,1:3,sum,na.rm=T)

    for(p in 1:np){
      Cp<-array(Cpred[,p,],c(nsim,nf))/apply(Cpred[,p,],1,sum)

      if(prod(round(CatchFrac[[p]],3)/round(Cp,3))!=1){
        print(Snames[p])
        print(cbind(CatchFrac[[p]],rep(NaN,nsim),round(Cp,3)))
        warning("Possible problem in catch fraction calculations")
      }

    }
  }

  # # apply hyperstability / hyperdepletion
  # Check that depletion is correct
  if (checks) {
    if (prod(round(Depletion,2)/ round(D_specified,2)) != 1) {
      print(cbind(round(Depletion,4),rep(NaN,nsim), round(D_specified,4)))
      warning("Possible problem in depletion calculations")
    }
  }

  if(!silent) message("Calculating MSY reference points")  # Print a progress update

  for(p in 1:np){   # --- Calculate MSY references ----

    V<-apply(FMt[,p,,,],1:3,sum)
    V<-nlz(V,c(1,3),"max")
    MSYrefs <- sapply(1:nsim, DLMtool::optMSY_eq, StockPars[[p]]$M_ageArray, StockPars[[p]]$Wt_age, StockPars[[p]]$Mat_age,
                      V=V, maxage=maxage,
                      R0=StockPars[[p]]$R0, SRrel=StockPars[[p]]$SRrel, hs=StockPars[[p]]$hs, yr.ind=(nyears-1):nyears)

    StockPars[[p]]$MSY <- MSYrefs[1, ]  # record the MSY results (Vulnerable)
    StockPars[[p]]$FMSY <- MSYrefs[2, ]  # instantaneous FMSY (Vulnerable)
    StockPars[[p]]$SSBMSY <- MSYrefs[3, ]  # Spawning Stock Biomass at MSY
    StockPars[[p]]$SSBMSY_SSB0 <- MSYrefs[4, ] # SSBMSY relative to unfished (SSB)
    StockPars[[p]]$BMSY_B0 <- MSYrefs[5, ] # Biomass relative to unfished (B0)
    StockPars[[p]]$BMSY <- MSYrefs[6,] # total biomass at MSY
    StockPars[[p]]$VBMSY <- (MSYrefs[1, ]/(1 - exp(-MSYrefs[2, ])))  # Biomass at MSY (Vulnerable)
    StockPars[[p]]$UMSY <- MSYrefs[1, ]/StockPars[[p]]$VBMSY  # exploitation rate [equivalent to 1-exp(-FMSY)]
    StockPars[[p]]$FMSY_M <- StockPars[[p]]$FMSY/StockPars[[p]]$M  # ratio of true FMSY to natural mortality rate M

  } # end of stocks

  # --- Code for deriving low biomass ----
  # (SSB where it takes MGThorizon x MGT to reach Bfrac of BMSY)
  # Znow<-apply(Z[,,nyears,]*N[,,nyears,],1:2,sum)/apply(N[,,nyears,],1:2,sum)
  # MGTsurv<-t(exp(-apply(Znow,1,cumsum)))
  # MGT<-apply(Agearray*(Mat_age[,,nyears]*MGTsurv),1,sum)/apply(Mat_age[,,nyears]*MGTsurv,1,sum)

  if(CalcBlow){
    #if(!silent) message("Calculating B-low reference points")              # Print a progress update

    #MGThorizon<-floor(HZN*MGT)

    #Blow <- sapply(1:nsim,getBlow, N, Asize, MSYrefs[3,],SSBpR, MPA, SSB0, nareas, retA, MGThorizon,
    #               Find,Perr_y,M_ageArray,hs,Mat_age, Wt_age,R0a,V,nyears,maxage,mov,
    #               Spat_targ,SRrel,aR,bR,Bfrac, maxF)
  }else{
    for(p in 1:np)StockPars[[p]]$Blow<-StockPars[[p]]$MGT<-rep(NA,nsim) # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

  }

  # --- Calculate Reference Yield ----
  if(!silent) message("Calculating reference yield - best fixed F strategy")  # Print a progress update

  for(p in 1:np)StockPars[[p]]$RefY <-StockPars[[p]]$MSY

  #  sapply(1:nsim, getFref3, Asize, nareas, maxage, N=N[,,nyears,, drop=FALSE], pyears=proyears,
  #               M_ageArray=M_ageArray[,,(nyears):(nyears+proyears)], Mat_age[,,(nyears):(nyears+proyears)],
  #               Wt_age=Wt_age[,,nyears:(nyears+proyears)],
  #               V=retA[, , (nyears + 1):(nyears + proyears), drop=FALSE],
  #               retA=retA[, , (nyears + 1):(nyears + proyears), drop=FALSE],
  #               Perr=Perr_y[,(nyears):(nyears+maxage+proyears-1)], mov, SRrel, Find,
  #               Spat_targ, hs, R0a, SSBpR, aR, bR, MPA=MPA, maxF=maxF, SSB0=SSB0)

  Ctemp<-array(NA,c(nsim,np,nf,maxage,nyears,nareas))
  CNind<-TEG(dim(Ctemp))
  Nind<-CNind[,c(1,2,4,5,6)]  # sim, stock, maxage, nyears, nareas

  # --- Calculate catch-at-age ----
  Ctemp[CNind]<-Biomass[Nind]*(1-exp(-Z[Nind]))*(FM[CNind]/Z[Nind])
  CB<-Ctemp# apply(Ctemp,1:5,sum)

  # --- Calculate retained-at-age ----

  Ctemp[CNind]<-N[Nind]*(1-exp(-Z[Nind]))*(FMret[CNind]/Z[Nind])
  Cret<-apply(Ctemp,1:5,sum)
  Cret[is.na(Cret)] <- 0

  Ctemp[CNind]<-Biomass[Nind]*(1-exp(-Z[Nind]))*(FMret[CNind]/Z[Nind])
  CBret<-Ctemp#apply(Ctemp,1:5,sum)

  # --- Calculate dead discarded-at-age ----
  CBdisc <- CB - CBret # discarded biomass

  # --- Simulate observed catch ----

  #if (length(control$Cbias_yr) ==0) {
    Cbiasa <- array(1, c(nsim,np,nf,nyears + proyears))  # Bias array
    for(p in 1:np){
      for(f in 1:nf){
        Cbiasa[,p,f,]<-ObsPars[[p]][[f]]$Cbias
      }
    }
  #} else {
  #  Cbiasa <- matrix(1, nsim, nyears+proyears)
  #  Cbiasa[,control$yrs] <- control$Cbias_yr
  #}

  Cerr <- array(NA,c(nsim,np,nf,nyears+proyears))

  for(p in 1:np){for(f in 1:nf){

    ObsPars[[p]][[f]]$Cerr<- Cerr[,p,f,] <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(ObsPars[[p]][[f]]$Csd, (nyears + proyears))),
                         sdconv(1, rep(ObsPars[[p]][[f]]$Csd, nyears + proyears))), c(nsim, nyears + proyears))
  }}

  # composite of bias and observation error
  for(p in 1:np){for(f in 1:nf){
      ObsPars[[p]][[f]]$Cobs <- Cbiasa[,p,f, 1:nyears] * Cerr[,p,f, 1:nyears] * apply(CBret[,p,f,,,], c(1,3), sum)  # Simulated observed retained catch (biomass)
  }}

  # --- Simulate observed catch-at-age ----

  # generate CAA from retained catch-at-age
  # CAA <- array(NA, dim = c(nsim, nf,np,nyears, maxage))  # Catch  at age array
  # cond <- apply(Cret, 1:2, sum, na.rm = T) < 1  # this is a fix for low sample sizes. If Cret is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
  # fixind <- as.matrix(cbind(expand.grid(1:nsim, 1:nyears), rep(floor(maxage/3), nyears)))  # more fix
  # Cret[fixind[cond, ]] <- 1  # puts a catch in the most vulnerable age class

  # a multinomial observation model for catch-at-age data
  for(p in 1:np){
    for(f in 1:nf){
      ObsPars[[p]][[f]]$CAA<-array(NA, dim = c(nsim, nyears, maxage))
      for (i in 1:nsim) {

        for (j in 1:nyears) {
          if (!sum( Cret[i, p,f,,j])) {
            ObsPars[[p]][[f]]$CAA[i, j, ] <- 0
          } else {
            ObsPars[[p]][[f]]$CAA[i, j, ] <- ceiling(-0.5 + rmultinom(1, ObsPars[[p]][[f]]$CAA_ESS[i], Cret[i, p,f,,j]) * ObsPars[[p]][[f]]$CAA_nsamp[i]/ObsPars[[p]][[f]]$CAA_ESS[i])
          }
        }
      }
    }
  }

  # --- Simulate observed catch-at-length ----
  # a multinomial observation model for catch-at-length data
  # assumed normally-distributed length-at-age truncated at 2 standard deviations from the mean

   for(p in 1:np){for(f in 1:nf){

    ObsPars[[p]][[f]]$CAL<- array(NA, c(nsim,  nyears, StockPars[[p]]$nCALbins))
    vn <- (apply(N[,p,,,], c(1,2,3), sum) * FleetPars[[p]][[f]]$retA[,,1:nyears]) # numbers at age that would be retained
    vn <- aperm(vn, c(1,3, 2))
    tempSize <- lapply(1:nsim, DLMtool::genSizeCompWrap, vn, StockPars[[p]]$CAL_binsmid, FleetPars[[p]][[f]]$retL, ObsPars[[p]][[f]]$CAL_ESS, ObsPars[[p]][[f]]$CAL_nsamp,
                       StockPars[[p]]$Linfarray,  StockPars[[p]]$Karray,  StockPars[[p]]$t0array, StockPars[[p]]$LenCV, truncSD=2.5)
    ObsPars[[p]][[f]]$CAL<-aperm(array(as.numeric(unlist(tempSize, use.names=FALSE)), dim=c(nyears, length(StockPars[[p]]$CAL_binsmid), nsim)), c(3,1,2))
    ObsPars[[p]][[f]]$LFC <- unlist(lapply(tempSize, function(x)  DLMtool::getfifth(x[nyears, ], StockPars[[p]]$CAL_binsmid)))
    ObsPars[[p]][[f]]$LFC[is.na(ObsPars[[p]][[f]]$LFC)] <- 1
    ObsPars[[p]][[f]]$LFC[ObsPars[[p]][[f]]$LFC<1] <- 1

  }}


  # calculate LFC

  # --- Simulate index of abundance from total biomass ----

  for(p in 1:np){for(f in 1:nf){
    ObsPars[[p]][[f]]$Ierr <- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(ObsPars[[p]][[f]]$Isd, nyears + proyears)),
                         sdconv(1, rep(ObsPars[[p]][[f]]$Isd, nyears + proyears))), c(nsim, nyears + proyears))
    II <- (apply(Biomass[,p,,,], c(1, 3), sum)^ObsPars[[p]][[f]]$betas) *  ObsPars[[p]][[f]]$Ierr[, 1:nyears]  # apply hyperstability / hyperdepletion
    ObsPars[[p]][[f]]$II <- II/apply(II, 1, mean)  # normalize
  }}

  for(p in 1:np){
    # --- Calculate vulnerable and spawning biomass abundance ----

    if (nsim > 1) A <- apply(VBiomass[,p, , nyears, ], 1, sum)  + apply(CB[,p,f , , nyears, ], 1, sum,na.rm=T) # Abundance before fishing
    if (nsim == 1) A <- sum(VBiomass[,p, , nyears, ]) +  sum(CB[,p,f,,nyears,],na.rm=T) # Abundance before fishing
    if (nsim > 1) Asp <- apply(SSB[,p, , nyears, ], 1, sum)  # SSB Abundance
    if (nsim == 1) Asp <- sum(SSB[,p, , nyears, ])  # SSB Abundance
    StockPars[[p]]$OFLreal <- A * StockPars[[p]]$FMSY  # the true simulated Over Fishing Limit
    StockPars[[p]]$A<-A
    StockPars[[p]]$Asp<-Asp

    # --- Simulate observed values in reference SBMSY/SB0 ----
    for(f in 1:nf){
      I3 <- apply(Biomass[,p,,,], c(1, 3), sum)^ObsPars[[p]][[f]]$betas  # apply hyperstability / hyperdepletion
      I3 <- I3/apply(I3, 1, mean)  # normalize index to mean 1
      # Iref <- apply(I3[, 1:5], 1, mean) * BMSY_B0  # return the real target abundance index corresponding to BMSY
      if (nsim > 1) Iref <- apply(I3[, 1:5], 1, mean) * StockPars[[p]]$SSBMSY_SSB0  # return the real target abundance index corresponding to BMSY
      if (nsim == 1) Iref <- mean(I3[1:5]) * StockPars[[p]]$SSBMSY_SSB0
      ObsPars[[p]][[f]]$I3<-I3
      ObsPars[[p]][[f]]$Iref<-Iref
    }
  }

  # --- Simulate observed values in steepness ----
  for(p in 1:np){for(f in 1:nf){

    if (!is.null(SampCpars[[p]][[f]][['hsim']])) {
      hsim <- SampCpars$hsim
      hbias <- hsim/hs  # back calculate the simulated bias
      if (Obs[[p]][[f]]@hbiascv == 0) hbias <- rep(1, nsim)
      ObsPars[[p]][[f]]$hbias <- hbias

    } else {
      if (is.null(SampCpars[[p]][[f]][['l_hbias']])) {
        hsim <- rep(NA, nsim)
        cond <- StockPars[[p]]$hs > 0.6
        hsim[cond] <- 0.2 + rbeta(sum(StockPars[[p]]$hs > 0.6), alphaconv((StockPars[[p]]$hs[cond] - 0.2)/0.8, (1 - (StockPars[[p]]$hs[cond] - 0.2)/0.8) * Obs[[p]][[f]]@hbiascv),
                                  betaconv((StockPars[[p]]$hs[cond] - 0.2)/0.8,  (1 - (StockPars[[p]]$hs[cond] - 0.2)/0.8) * Obs[[p]][[f]]@hbiascv)) * 0.8

        hsim[!cond] <- 0.2 + rbeta(sum(StockPars[[p]]$hs <= 0.6), alphaconv((StockPars[[p]]$hs[!cond] - 0.2)/0.8,  (StockPars[[p]]$hs[!cond] - 0.2)/0.8 * Obs[[p]][[f]]@hbiascv),
                                   betaconv((StockPars[[p]]$hs[!cond] - 0.2)/0.8, (StockPars[[p]]$hs[!cond] - 0.2)/0.8 * Obs[[p]][[f]]@hbiascv)) * 0.8

        hbias <- hsim/StockPars[[p]]$hs  # back calculate the simulated bias
        if (Obs[[p]][[f]]@hbiascv == 0) hbias <- rep(1, nsim)
        ObsPars[[p]][[f]]$hbias <- hbias
      } else {
        l_hbias <- sample(SampCpars[[p]][[f]]$l_hbias, nsim, replace=TRUE)

        P <- (StockPars[[p]]$hs-.2)/0.8
        hs_logit <- log(P/(1-P))
        P2 <- (exp(hs_logit)* l_hbias)/(1+exp(hs_logit) * l_hbias)
        P2[is.nan(P2)] <- 1
        hsim <- (P2 * 0.8) + 0.2
        hbias <- hsim/hs  # back calculate the simulated bias
        ObsPars[[p]][[f]]$hbias <- hbias

      }
    }  # end of is hsim in cpars
  }} # end of p f


  # --- Simulate error in observed recruitment index ----
  for(p in 1:np){for(f in 1:nf){
    ObsPars[[p]][[f]]$Recerr<- array(rlnorm((nyears + proyears) * nsim, mconv(1, rep(ObsPars[[p]][[f]]$Recsd, (nyears + proyears))),
                         sdconv(1, rep(ObsPars[[p]][[f]]$Recsd, nyears + proyears))), c(nsim, nyears + proyears))

    # --- Simulate observation error in BMSY/B0 ----
    ntest <- 20  # number of trials
    BMSY_B0bias <- array(rlnorm(nsim * ntest, mconv(1, Obs[[p]][[f]]@BMSY_B0biascv), sdconv(1, Obs[[p]][[f]]@BMSY_B0biascv)), dim = c(nsim, ntest))  # trial samples of BMSY relative to unfished
    # test <- array(BMSY_B0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0
    test <- array(StockPars[[p]]$SSBMSY_SSB0 * BMSY_B0bias, dim = c(nsim, ntest))  # the simulated observed BMSY_B0
    indy <- array(rep(1:ntest, each = nsim), c(nsim, ntest))  # index

    # indy[test > 0.9] <- NA  # interval censor
    indy[test > max(0.9, max(StockPars[[p]]$SSBMSY_SSB0))] <- NA  # interval censor

    BMSY_B0bias <- BMSY_B0bias[cbind(1:nsim, apply(indy, 1, min, na.rm = T))]  # sample such that BMSY_B0<90%
    ObsPars[[p]][[f]]$BMSY_B0bias <- BMSY_B0bias
  }}
  # --- Implementation error time series ----

  for(p in 1:np){for(f in 1:nf){
      ImpPars[[p]][[f]]$TAC_f <- array(rlnorm(proyears * nsim, mconv(ImpPars[[p]][[f]]$TACFrac, ImpPars[[p]][[f]]$TACSD),
                            sdconv(ImpPars[[p]][[f]]$TACFrac, ImpPars[[p]][[f]]$TACSD)), c(nsim, proyears))  # composite of TAC fraction and error

      ImpPars[[p]][[f]]$E_f <- array(rlnorm(proyears * nsim, mconv(ImpPars[[p]][[f]]$TAEFrac, ImpPars[[p]][[f]]$TAESD),
                          sdconv(ImpPars[[p]][[f]]$TAEFrac, ImpPars[[p]][[f]]$TAESD)), c(nsim, proyears))  # composite of TAC fraction and error

      ImpPars[[p]][[f]]$SizeLim_f<-array(rlnorm(proyears * nsim, mconv(ImpPars[[p]][[f]]$SizeLimFrac, ImpPars[[p]][[f]]$SizeLimSD),
                              sdconv(ImpPars[[p]][[f]]$SizeLimFrac, ImpPars[[p]][[f]]$SizeLimSD)), c(nsim, proyears))  # composite of TAC fraction and error
  }}

  # --- Populate Data object with Historical Data ----
  DataList<-new('list')

  for(p in 1:np)DataList[[p]]<-new('list')

  for(p in 1:np){for(f in 1:nf){

    Data <- new("Data", stock = "MSE")  # create a blank DLM data object
    if (MOM@reps == 1) Data <- DLMtool::OneRep(Data)  # make stochastic variables certain for only one rep
    Data <- replic8(Data, nsim)  # make nsim sized slots in the DLM data object
    Data@Name <- paste(Stocks[[p]]@Name,Fleets[[p]][[f]]@Name,sep="-")
    Data@Year <- 1:nyears
    Data@Cat <-  ObsPars[[p]][[f]]$Cobs
    Data@Ind <-  ObsPars[[p]][[f]]$II
    Data@Rec <- apply(N[, p,1, , ], c(1, 2), sum) *  ObsPars[[p]][[f]]$Recerr[, 1:nyears]
    Data@t <- rep(nyears, nsim)
    Data@AvC <- apply(ObsPars[[p]][[f]]$Cobs, 1, mean)
    Data@Dt <- ObsPars[[p]][[f]]$Dbias * StockPars[[p]]$Depletion * rlnorm(nsim, mconv(1,  ObsPars[[p]][[f]]$Derr), sdconv(1,  ObsPars[[p]][[f]]$Derr))
    Data@Mort <- StockPars[[p]]$M * ObsPars[[p]][[f]]$Mbias
    Data@FMSY_M <- StockPars[[p]]$FMSY_M * ObsPars[[p]][[f]]$FMSY_Mbias
    # Data@BMSY_B0 <- BMSY_B0 * BMSY_B0bias
    Data@BMSY_B0 <- StockPars[[p]]$SSBMSY_SSB0 * ObsPars[[p]][[f]]$BMSY_B0bias
    Data@Cref <- StockPars[[p]]$MSY * ObsPars[[p]][[f]]$Crefbias
    Data@Bref <- StockPars[[p]]$VBMSY * ObsPars[[p]][[f]]$Brefbias
    Data@Iref <- ObsPars[[p]][[f]]$Iref * ObsPars[[p]][[f]]$Irefbias
    Data@LFC <- ObsPars[[p]][[f]]$LFC * ObsPars[[p]][[f]]$LFCbias
    Data@LFS <- FleetPars[[p]][[f]]$LFS[nyears,] * ObsPars[[p]][[f]]$LFSbias
    Data@CAA <- ObsPars[[p]][[f]]$CAA
    Data@Dep <- ObsPars[[p]][[f]]$Dbias * StockPars[[p]]$Depletion * rlnorm(nsim, mconv(1,  ObsPars[[p]][[f]]$Derr), sdconv(1,  ObsPars[[p]][[f]]$Derr))
    Data@Abun <- StockPars[[p]]$A * ObsPars[[p]][[f]]$Abias * rlnorm(nsim, mconv(1, ObsPars[[p]][[f]]$Aerr), sdconv(1, ObsPars[[p]][[f]]$Aerr))
    Data@SpAbun <- StockPars[[p]]$Asp * ObsPars[[p]][[f]]$Abias * rlnorm(nsim, mconv(1, ObsPars[[p]][[f]]$Aerr), sdconv(1, ObsPars[[p]][[f]]$Aerr))
    Data@vbK <- StockPars[[p]]$K * ObsPars[[p]][[f]]$Kbias
    Data@vbt0 <- StockPars[[p]]$t0 * ObsPars[[p]][[f]]$t0bias
    Data@LenCV <- StockPars[[p]]$LenCV # * LenCVbias
    Data@vbLinf <- StockPars[[p]]$Linf * ObsPars[[p]][[f]]$Linfbias
    Data@L50 <- StockPars[[p]]$L50 * ObsPars[[p]][[f]]$lenMbias
    Data@L95 <- StockPars[[p]]$L95 * ObsPars[[p]][[f]]$lenMbias
    Data@L95[Data@L95 > 0.9 * Data@vbLinf] <- 0.9 * Data@vbLinf[Data@L95 > 0.9 * Data@vbLinf]  # Set a hard limit on ratio of L95 to Linf
    Data@L50[Data@L50 > 0.9 * Data@L95] <- 0.9 * Data@L95[Data@L50 > 0.9 * Data@L95]  # Set a hard limit on ratio of L95 to Linf
    Data@steep <- StockPars[[p]]$hs * ObsPars[[p]][[f]]$hbias
    Data@sigmaR <- StockPars[[p]]$procsd # sigmaR assumed no obs error
    Data@CAL_bins <- StockPars[[p]]$CAL_bins
    Data@CAL <- ObsPars[[p]][[f]]$CAL
    MLbin <- (StockPars[[p]]$CAL_bins[1:(length(StockPars[[p]]$CAL_bins) - 1)] + StockPars[[p]]$CAL_bins[2:length(StockPars[[p]]$CAL_bins)])/2
    temp <- ObsPars[[p]][[f]]$CAL * rep(MLbin, each = nsim * nyears)
    Data@ML <- apply(temp, 1:2, sum)/apply(ObsPars[[p]][[f]]$CAL, 1:2, sum)
    Data@Lc <- array(MLbin[apply(ObsPars[[p]][[f]]$CAL, 1:2, which.max)], dim = c(nsim, nyears))
    nuCAL <- ObsPars[[p]][[f]]$CAL
    for (i in 1:nsim) for (j in 1:nyears) nuCAL[i, j, 1:match(max(1, Data@Lc[i, j]), MLbin, nomatch=1)] <- NA
    temp <- nuCAL * rep(MLbin, each = nsim * nyears)
    Data@Lbar <- apply(temp, 1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE)
    Data@MaxAge <- maxage
    Data@Units <- "unitless"
    Data@Ref <- StockPars[[p]]$OFLreal
    Data@Ref_type <- "Simulated OFL"
    Data@wla <- rep(StockPars[[p]]$a, nsim)
    Data@wlb <- rep(StockPars[[p]]$b, nsim)
    Data@nareas <- nareas

    # put all the operating model parameters in one table

    OMtable<-cbind(
      as.data.frame(StockPars[[p]][c('RefY',"M","D","A","SSBMSY_SSB0","FMSY_M","Msd","procsd","MSY",
                                   "FMSY",'Linf','K','t0','hs','Linfsd','Ksd','OFLreal','Size_area_1',
                                   'Frac_area_1','Prob_staying','AC','L50','L95','B0','N0','SSB0','BMSY_B0','Blow','MGT','BMSY','SSBMSY','Mexp','Fdisc')]),

      as.data.frame(FleetPars[[p]][[f]][c('Esd','dFfinal','qinc','qcv','Spat_targ')]),

      as.data.frame(ImpPars[[p]][[f]][c('TACSD','TACFrac','TAESD','TAEFrac','SizeLimSD','SizeLimFrac')]),

      data.frame(ageM=StockPars[[p]]$ageM[,nyears], L5=FleetPars[[p]][[f]]$L5[nyears,],LFS=FleetPars[[p]][[f]]$LFS[nyears,],
               Vmaxlen=FleetPars[[p]][[f]]$Vmaxlen[nyears,], LFC=ObsPars[[p]][[f]]$LFC,
               LR5=FleetPars[[p]][[f]]$LR5[nyears,],  LFR=FleetPars[[p]][[f]]$LFR[nyears,],
               Rmaxlen=FleetPars[[p]][[f]]$Rmaxlen[nyears,], DR=FleetPars[[p]][[f]]$DR[nyears,])

    )

    OMtable <- OMtable[,order(names(OMtable))]
    Data@OM <- OMtable

    ObsTable <- as.data.frame(ObsPars[[p]][[f]][1:24])
    ObsTable <- ObsTable[,order(names(ObsTable))]
    Data@Obs <- ObsTable # put all the observation error model parameters in one table

    Data@LHYear <- nyears  # Last historical year is nyears (for fixed MPs)
    histCatches <- apply(CBret[,p,f,,,], c(1, 3), sum,na.rm=T)
    Data@MPrec <- histCatches[, nyears]
    Data@MPeff <- rep(1, nsim)
    Data@Misc <- vector("list", nsim)
    DataList[[p]][[f]]<-Data

  }} # end of np nf

  # --- Return Historical Simulations and Data from last historical year ----
  if (Hist) { # Stop the model after historical simulations are complete
    if(!silent) message("Returning historical simulations")

    nout <- apply(N, c(1,2, 4), sum)
    vb <- apply(VBiomass, c(1,2,4), sum)
    b <- apply(Biomass, c(1,2,4), sum)
    ssb <- apply(SSB, c(1,2,4), sum)
    Cc <- apply(CB, c(1,2,3,5), sum)
    rec <- apply((N)[,,1, , ], c(1,2,3), sum)

    TSdata <- list(VB=vb, SSB=ssb, Bio=b, Catch=Cc, Rec=rec, N=nout)
    AtAge <- list(Z=Z, FM=FM, FMret=FMret)

    StockPars$Depletion <- Depletion
    FleetPars$qs <- qs
    SampPars <- c(StockPars, FleetPars, ObsPars, ImpPars)
    Data@Misc <- list()
    HistData <- list(StockPars=StockPars, FleetPars=FleetPars, ObsPars=ObsPars, ImpPars=ImpPars,
                     TSdata=TSdata, AtAge=AtAge, Stocks=Stocks, Fleets=Fleets, Obs=Obs, Imps=Imps,
                     Data=DataList)
    return(HistData)
  }

  # assign('Data',Data,envir=.GlobalEnv) # for debugging fun

  # Detecting MP specification  -----------------------------------------------------------------------------------------------------------------------------------

  if(identical(ldim(MPs),ldim(Fleets))){
    message("Byfleet mode: you have specified an MP for each stock and fleet. Only fleet-specific data (e.g. catches and indices) will be used to set advice for each fleet for each stock")
    MPcond<-"byfleet"
    nMP <- length(MPs[[1]][[1]])
    MPrefs<-array(NA,c(nMP,nf,np))
    MPrefs[]<-unlist(MPs)
  }else if(np==1&nf==1){
    nMP <- length(MPs[[1]][[1]])
    MPcond<-"bystock"
    message("runMSE checking: you have specified a single stock and fleet. For analysis you should be using runMSE(). Use this only for debugging against runMSE.")

    MPrefs<-array(NA,c(nMP,nf,np))
    MPrefs[]<-unlist(MPs)
  }else{
    if(ldim(MPs)==ldim(Fleets)[1]){ # not a two-tier list
      message("Bystock mode: you have specified a vector of MPs for each stock, but not a vector of MPs for each stock and fleet. The catch data for these fleets will be combined, a single MP will be used to set a single TAC for all fleets combined that will be allocated between the fleets according to recent catches")
      MPcond<-"bystock"
      nMP<-length(MPs[[1]])
      MPrefs<-array(NA,c(nMP,nf,np))
      for(p in 1:np)MPrefs[,,p]<-MPs[[p]]
    }
    if(class(MPs)!="list"){
      if(class(get(MPs[1]))=="MMP"){
        message("MMP mode: you have specified multi-fleet, multi-stock MPs of class MMP. This class of MP accepts all data objects (stocks x fleets) to simultaneously make a recommendation specific to each stock and fleet")
        MPcond<-"MMP"
        nMP<-length(MPs)
        MPrefs<-array(NA,c(nMP,nf,np))
        MPrefs[]<-MPs
      }else if(class(get(MPs[1]))=="MP"){
        message("Complex mode: you have specified a vector of MPs rather than a list of MPs, one list position for MP type. The same MP will be applied to the aggregate data for all stocks and fleets. The MP will, for example, be used to set a single TAC for all stocks and fleets combined. This will be allocated among fleets according to recent catches and among stocks according to available, vulnerable biomass")
        MPcond<-"complex"
        MPtemp<-MPs
        nMP<-length(MPs)
        MPrefs<-array(NA,c(nMP,nf,np))
        MPrefs[]<-unlist(MPs)
      }
    } # not a list
  } # end of two

  if(class(MPs)=="list"){
    allMPs<-unlist(MPs)
  }else{
    allMPs<-MPs
  }
    # --- Check MPs ----
  if (CheckMPs & MPcond != "MMP") {
    if(!silent) message("Determining available methods")  # print an progress report
    PosMPs <- Can( Data, timelimit = timelimit)  # list all the methods that could be applied to a DLM data object

    if (!is.na(allMPs[1])) {
      cant <- allMPs[!allMPs %in% PosMPs]
      if (length(cant) > 0) {
        if(!silent) stop(paste0("Cannot run some MPs:", DLMtool::DLMdiag(Data, "not available", funcs1=cant, timelimit = timelimit)))
      }
    }
  }

  # Create a data object for each method (they have identical historical data and branch in projected years)
  # also create the CALout (true catch at length) by stock and fleet as a list since nCALbins varies among stocks (ragged)
  MSElist<-CALout<-list('list')
  for(p in 1:np){
    MSElist[[p]]<-new('list')
    CALout[[p]]<-new('list')
    for(f in 1:nf){
      MSElist[[p]][[f]]<-list(DataList[[p]][[f]])[rep(1, nMP)]
      CALout[[p]][[f]]<-list()
  }}

  B_BMSYa <- Ba <- SSBa <- VBa <- array(NA, dim = c(nsim,np,nMP, proyears))  # store the projected B_BMSY
  FMa <-F_FMSYa<- Ca <- CaRet <- TACa <- Effort <- array(NA, dim = c(nsim,np,nf, nMP, proyears))  # store the projected fishing mortality rate
  CAAout <- array(NA, dim = c(nsim,np,nf, nMP, maxage))  # store the population-at-age in last projection year
  PAAout <-  array(NA, dim = c(nsim,np,nMP,maxage))
  #CALout <- array(NA, dim = c(nsim,np,nf, nMP, nCALbins))  # store the population-at-length in last projection year
  # SPRa <- array(NA,dim=c(nsim,nMP,proyears)) # store the Spawning Potential Ratio

  # --- Calculate MSY statistics for each projection year ----
  MSY_P <- FMSY_P <- SSBMSY_P <- array(NA, dim=c(nsim,np, nMP, proyears))

  for(p in 1:np){
    MSY_P[,p,,] <- StockPars[[p]]$MSY
    FMSY_P[,p,,] <- StockPars[[p]]$FMSY
    SSBMSY_P[,p,,] <- StockPars[[p]]$SSBMSY
  }

  interval<-MOM@interval
  if (length(MOM@interval) != nMP) interval <- rep(interval, nMP)[1:nMP]
  if (!all(interval == interval[1])) {
    message("Variable management intervals:")
    df <- data.frame(MP=MPs,interval=interval)
    message(paste(capture.output(print(df)), collapse = "\n"))
  }

  # --- Begin loop over MPs ----
  mm <- 1 # for debugging
  TAC_A<-array(NA,c(nsim,np,nf)) # Temporary store of the TAC
  TAE_A<-array(NA,c(nsim,np,nf)) # Temporary store of the TAE
  MPrecs_A_blank<-list() # Temporary Hierarcical list of MPrec objects
  for(p in 1:np)MPrecs_A_blank[[p]]<-list()
  LastTAE<- histTAE<- Effort_pot<-LastAllocat<-LastCatch<-TACused<-array(NA,c(nsim,np,nf))
  LastSpatial<-array(NA,c(nareas,np,nf,nsim))
  V_Pt<-array(NA,c(nsim,nf,maxage,nyears+proyears)) # temporary vulnerability for MSY calcs combined over fleets

  # ===========================================================================================================================================================================
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # ===========================================================================================================================================================================

  mm<-1

  for (mm in 1:nMP) {  # MSE Loop over methods

    if(!silent){
      message(" ----- ", mm, "/", nMP, " MPs, Running MSE for: ")  # print a progress report
      for(p in 1:np){
        MPrep<-data.frame(MPrefs[mm,,p])
        row.names(MPrep)<-Fnames[,p]
        names(MPrep)=Snames[p]
        print(MPrep)
      }
      message(" --------------------------------- ")
    }

    checkNA <- array(0,c(np,nf,proyears)) # save number of NAs

    for(p in 1:np){
      for(f in 1:nf){

        # reset selectivity parameters for projections
        FleetPars[[p]][[f]]$L5_P <- FleetPars[[p]][[f]]$L5
        FleetPars[[p]][[f]]$LFS_P <- FleetPars[[p]][[f]]$LFS
        FleetPars[[p]][[f]]$Vmaxlen_P <- FleetPars[[p]][[f]]$Vmaxlen
        FleetPars[[p]][[f]]$SLarray_P <- FleetPars[[p]][[f]]$SLarray # selectivity at length array - projections
        FleetPars[[p]][[f]]$V_P <- FleetPars[[p]][[f]]$V  #  selectivity at age array - projections

        # reset retention parametersfor projections
        FleetPars[[p]][[f]]$LR5_P <- FleetPars[[p]][[f]]$LR5
        FleetPars[[p]][[f]]$LFR_P <- FleetPars[[p]][[f]]$LFR
        FleetPars[[p]][[f]]$Rmaxlen_P <- FleetPars[[p]][[f]]$Rmaxlen
        FleetPars[[p]][[f]]$retA_P <- FleetPars[[p]][[f]]$retA # retention at age array - projections
        FleetPars[[p]][[f]]$retL_P <- FleetPars[[p]][[f]]$retL # retention at length array - projections
        FleetPars[[p]][[f]]$DR_P <- FleetPars[[p]][[f]]$DR # Discard ratio for projections

        FleetPars[[p]][[f]]$FM_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
        FleetPars[[p]][[f]]$FM_Pret <- array(NA, dim = c(nsim, maxage, proyears, nareas)) # retained F
        FleetPars[[p]][[f]]$FM_nospace <- array(NA, dim = c(nsim, maxage, proyears, nareas))  # stores prospective F before reallocation to new areas
        FleetPars[[p]][[f]]$FML <- array(NA, dim = c(nsim, nareas))  # last apical F
        FleetPars[[p]][[f]]$Z_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
        FleetPars[[p]][[f]]$CB_P <- array(NA, dim = c(nsim,maxage, proyears, nareas))
        FleetPars[[p]][[f]]$CB_Pret <- array(NA, dim = c(nsim,maxage, proyears, nareas)) # retained catch

      }

      StockPars[[p]]$Fdisc_P <- StockPars[[p]]$Fdisc # Discard mortality for projectons
      StockPars[[p]]$N_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
      StockPars[[p]]$Biomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
      StockPars[[p]]$VBiomass_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))
      StockPars[[p]]$SSN_P <-array(NA, dim = c(nsim,maxage, proyears, nareas))
      StockPars[[p]]$SSB_P <- array(NA, dim = c(nsim, maxage, proyears, nareas))

    }

    N_P <- array(NA, dim = c(nsim, np, maxage, proyears, nareas))
    Biomass_P <- array(NA, dim = c(nsim,np, maxage, proyears, nareas))
    VBiomass_P <- array(NA, dim = c(nsim,np, maxage, proyears, nareas))
    SSN_P <-array(NA, dim = c(nsim,np, maxage, proyears, nareas))
    SSB_P <- array(NA, dim = c(nsim,np, maxage, proyears, nareas))
    FMt_P <- array(NA, dim = c(nsim, np, maxage, proyears, nareas))
    Z_P <- array(NA, dim = c(nsim, np, maxage, proyears, nareas))
    FM_P <- array(NA, dim = c(nsim, np,nf,maxage, proyears, nareas))
    FMret_P <- array(NA, dim = c(nsim,np,nf, maxage, proyears, nareas))
    VBF_P<-array(NA, dim = c(nsim,np,nf, maxage, proyears, nareas))

    # indexes
    SAYRL <- as.matrix(expand.grid(1:nsim, 1:maxage, nyears, 1:nareas))  # Final historical year
    SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, 1 + nyears, 1:nareas))  # Trajectory year
    SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, 1, 1:nareas))
    SYt <- SAYRt[, c(1, 3)]
    SAYt <- SAYRt[, 1:3]
    SR <- SAYR[, c(1, 4)]
    SA1 <- SAYR[, 1:2]
    S1 <- SAYR[, 1]
    SY1 <- SAYR[, c(1, 3)]
    SAY1 <- SAYRt[, 1:3]
    SYA <- as.matrix(expand.grid(1:nsim, 1, 1:maxage))  # Projection year
    SY <- SYA[, 1:2]
    SA <- SYA[, c(1, 3)]
    SAY <- SYA[, c(1, 3, 2)]
    S <- SYA[, 1]

    # -- First projection year ----
    y <- 1
    Perr<-hs<-R0<-SRrel<-K<-Linf<-t0<-M<-array(NA,c(nsim,np))
    aR<-bR<-R0a<-SSBpR<-Asize<-array(NA,c(nsim,np,nareas))
    mov<-array(NA,c(nsim,np,maxage,nareas,nareas,nyears+proyears))
    Spat_targ_y<-array(NA,c(nsim,np,nf))
    M_agecur_y<-Mat_agecur_y<-array(NA,c(nsim,np,maxage))
    a_y<-b_y<-rep(NA,np)

    for(p in 1:np){
      Perr[,p]<-StockPars[[p]]$Perr_y[,nyears+maxage-1]
      hs[,p]<-StockPars[[p]]$hs
      aR[,p,]<-StockPars[[p]]$aR
      bR[,p,]<-StockPars[[p]]$bR
      mov[,p,,,,]<-StockPars[[p]]$mov
      for(f in 1:nf)Spat_targ_y[,p,f]<-FleetPars[[p]][[f]]$Spat_targ
      SRrel[,p]<-StockPars[[p]]$SRrel
      M_agecur_y[,p,]<-StockPars[[p]]$M_ageArray[,,nyears]
      Mat_agecur_y[,p,]<-StockPars[[p]]$Mat_age[,,nyears]
      K[,p]<-StockPars[[p]]$Karray[,nyears]
      Linf[,p]<-StockPars[[p]]$Linfarray[,nyears]
      t0[,p]<-StockPars[[p]]$t0array[,nyears]
      M[,p]<-StockPars[[p]]$M
      R0[,p]<-StockPars[[p]]$R0
      R0a[,p,]<-StockPars[[p]]$R0a
      SSBpR[,p,]<-StockPars[[p]]$SSBpR
      a_y[p]<-StockPars[[p]]$a
      b_y[p]<-StockPars[[p]]$b
      Asize[,p,]<-StockPars[[p]]$Asize
    }

    NextYrN <- sapply(1:nsim, function(x)
     popdynOneMICE(np,nf,nareas, maxage, Ncur=array(N[x,,,nyears,],c(np,maxage,nareas)), Vcur=array(VF[x,,,,nyears],c(np,nf,maxage)),
                  #Retcur=array(FretA[x,,,,nyears],c(np,nf,maxage)),
                  #Fcur=apply(apply(array(FM[x,,,,nyears,],c(np,nf,maxage,nareas)),c(1,2,4),max),c(1,2),sum),   # note that Fcur is apical F but, in popdynOneMICE it is DIVIDED in future years between the two areas depending on vulnerabile biomass. So to get Fcur you need to sum over areas (a bit weird)
                  FMretx=array(FMret[x,,,,nyears,],c(np,nf,maxage,nareas)),
                  FMx=array(FM[x,,,,nyears,],c(np,nf,maxage,nareas)),   # note that Fcur is apical F but, in popdynOneMICE it is DIVIDED in future years between the two areas depending on vulnerabile biomass. So to get Fcur you need to sum over areas (a bit weird)
                  PerrYrp=Perr[x,], hsx=hs[x,], aRx=matrix(aR[x,,],nrow=np), bRx=matrix(bR[x,,],nrow=np),
                  movy=array(mov[x,,,,,nyears],c(np,maxage,nareas,nareas)), Spat_targ=array(Spat_targ_y[x,,],c(np,nf)), SRrelx=SRrel[x,],
                  M_agecur=matrix(M_agecur_y[x,,],nrow=np), Mat_agecur=matrix(Mat_agecur_y[x,,],nrow=np),
                  Asizex=matrix(Asize[x,,],ncol=nareas),Kx =K[x,], Linfx=Linf[x,],t0x=t0[x,],Mx=M[x,],
                  R0x=R0[x,],R0ax=matrix(R0a[x,,],nrow=np),SSBpRx=matrix(SSBpR[x,,],nrow=np),ax=a_y,
                  bx=b_y,Rel=Rel,SexPars=SexPars,x=x))
    # x=1; Ncur=array(N[x,,,nyears,],c(np,maxage,nareas)); Vcur=array(VF[x,,,,nyears],c(np,nf,maxage)); FMretx=array(FMret[x,,,,nyears,],c(np,nf,maxage,nareas));  FMx=array(FM[x,,,,nyears,],c(np,nf,maxage,nareas));   PerrYrp=Perr[x,]; hsx=hs[x,]; aRx=matrix(aR[x,,],nrow=np); bRx=matrix(bR[x,,],nrow=np);  movy=array(mov[x,,,,,nyears],c(np,maxage,nareas,nareas)); Spat_targ=array(Spat_targ_y[x,,],c(np,nf)); SRrelx=SRrel[x,]; M_agecur=matrix(M_agecur_y[x,,],nrow=np); Mat_agecur=matrix(Mat_agecur_y[x,,],nrow=np);   Asizex=matrix(Asize[x,,],ncol=nareas); Kx =K[x,]; Linfx=Linf[x,];t0x=t0[x,];Mx=M[x,]; R0x=R0[x,]; R0ax=matrix(R0a[x,,],nrow=np); SSBpRx=matrix(SSBpR[x,,],nrow=np); ax=a_y; bx=b_y; Rel=Rel; SexPars=SexPars; x=x

    N_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[1,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
    Biomass_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[23,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
    SSN_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[24,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
    SSB_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[25,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
    VBiomass_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[19,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
    FML <- apply(array(FM[, ,,, nyears, ],c(nsim,np,nf,maxage,nareas)), c(1, 3), max)

    #FM_P[,,,,1,] <- aperm(array(as.numeric(unlist(NextYrN[17,], use.names=FALSE)), dim=c(np ,nf,maxage, nareas, nsim)), c(5,1,2,3,4))
    #FMret_P[,,,,1,] <- aperm(array(as.numeric(unlist(NextYrN[18,], use.names=FALSE)), dim=c(np ,nf,maxage, nareas, nsim)), c(5,1,2,3,4))
    #VBF[,,,,1,] <- aperm(array(as.numeric(unlist(NextYrN[20,], use.names=FALSE)), dim=c(np ,nf,maxage, nareas, nsim)), c(5,1,2,3,4))
    #Z_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[21,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
    #FMt_P[,,,1,]<-aperm(array(as.numeric(unlist(NextYrN[22,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))

    for(p in 1:np){
      StockPars[[p]]$N_P<-N_P[,p,,,]
      StockPars[[p]]$Biomass_P<-Biomass_P[,p,,,]
      StockPars[[p]]$SSN_P<-SSN_P[,p,,,]
      StockPars[[p]]$SSB_P<-SSB_P[,p,,,]
      StockPars[[p]]$VBiomass_P<-VBiomass_P[,p,,,]
      #for(f in 1:nf)FleetPars[[p]][[f]]$FML<-FML[]
    }

    # -- apply MP in initial projection year ----
    # - Combined MP -

    if(MPcond=="MMP"){

      DataList<-getDataList(MSElist,mm) # returns a hierarchical list object stock then fleet of Data objects
      MPRecs_A <- applyMMP(DataList, MP = MPs[mm], reps = 1, silent=TRUE)  # # returns a hierarchical list object stock then fleet then slot type of Rec
      Data_p_A <- MPrecs_A_blank
      for(p in 1:np)for(f in 1:nf){
        Data_p_A[[p]][[f]]<-MSElist[[p]][[f]][[mm]]
        Data_p_A[[p]][[f]]@TAC<-MPRecs_A[[p]][[f]]$TAC # record TAC rec in Data
      }

    }else if(MPcond=="complex"){

      MPRecs_A <- Data_p_A <- MPrecs_A_blank # A temporary blank hierarchical list object stock by fleet
      realVB<-apply(VBiomass[,,,1:nyears,],c(1,2,4),sum,na.rm=T) # need this for aggregating data and distributing TACs over stocks

      curdat<-multiDataS(MSElist,StockPars,np,mm,nf,realVB)
      runMP <- applyMP(curdat, MPs = MPs[mm], reps = 1, silent=TRUE)  # Apply MP

      Stock_Alloc<-realVB[,,nyears]/apply(realVB[,,nyears],1,sum)

      for(p in 1:np)  for(f in 1:nf){
        MPRecs_A[[p]][[f]]<-runMP[[1]][[1]]
        MPRecs_A[[p]][[f]]$TAC<-runMP[[1]][[1]]$TAC*MOM@Allocation[[p]][,f]*Stock_Alloc[,p]
        MPRecs_A[[p]][[f]]$Effort<-runMP[[1]][[1]]$Effort*MOM@Efactor[[p]][,f]

        if(length(MPRecs_A[[p]][[f]]$Effort)>0) if(is.na(MPRecs_A[[p]][[f]]$Effort[1,1])) MPRecs_A[[p]][[f]]$Effort<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$Effort))
        if(length(MPRecs_A[[p]][[f]]$TAC)>0) if(is.na(MPRecs_A[[p]][[f]]$TAC[1,1])) MPRecs_A[[p]][[f]]$TAC<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))
        if(is.na(MPRecs_A[[p]][[f]]$Spatial[1,1])) MPRecs_A[[p]][[f]]$Spatial<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))

        Data_p_A[[p]][[f]]<-runMP[[2]]
        Data_p_A[[p]][[f]]@TAC<-MPRecs_A[[p]][[f]]$TAC
      }

    }else{

      MPRecs_A <- Data_p_A <- MPrecs_A_blank # A temporary blank hierarchical list object stock by fleet

      for(p in 1:np){

        if(MPcond=="bystock"){

          if(nf>1){
            curdat<-multiData(MSElist,StockPars,p,mm,nf)
          }else{
            curdat<-MSElist[[p]][[f]][[mm]]
          }

          runMP <- applyMP(curdat, MPs = MPs[[p]][mm], reps = 1, silent=TRUE)  # Apply MP

          # Do allocation calcs
          TAC_A[,p,]<-array(as.vector(unlist(runMP[[1]][[1]]$TAC))*MOM@Allocation[[p]],c(nsim,nf))
          TAE_A[,p,]<-array(as.vector(unlist(runMP[[1]][[1]]$Effort))*MOM@Efactor[[p]],c(nsim,nf))

          for(f in 1:nf){
            MPRecs_A[[p]][[f]]<-runMP[[1]][[1]]
            MPRecs_A[[p]][[f]]$TAC<-matrix(TAC_A[,p,f],nrow=1) # copy allocated TAC
            MPRecs_A[[p]][[f]]$Effort<-matrix(TAE_A[,p,f],nrow=1)
            # This next line is to make the NULL effort recommendations of an output control MP compatible with CalcMPdynamics (expects a null matrix)
            if(is.na(MPRecs_A[[p]][[f]]$Effort[1,1])) MPRecs_A[[p]][[f]]$Effort<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$Effort))
            if(is.na(MPRecs_A[[p]][[f]]$TAC[1,1])) MPRecs_A[[p]][[f]]$TAC<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))
            if(is.na(MPRecs_A[[p]][[f]]$Spatial[1,1])) MPRecs_A[[p]][[f]]$Spatial<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))

            Data_p_A[[p]][[f]]<-runMP[[2]]
            Data_p_A[[p]][[f]]@TAC<-MPRecs_A[[p]][[f]]$TAC   # copy allocated tAC
          }

        }else if(MPcond=="byfleet"){

          for(f in 1:nf){

            curdat<-MSElist[[p]][[f]][[mm]]
            runMP <- applyMP(curdat, MPs = MPrefs[p,f,mm], reps = 1, silent=TRUE)  # Apply MP
            MPRecs_A[[p]][[f]]<-runMP[[1]][[1]]
            Data_p_A[[p]][[f]]<-runMP[[2]]
            Data_p_A[[p]][[f]]@TAC <- MPRecs_A[[p]][[f]]$TAC

          }

        }

      } # end of stocks

    }

    for(p in 1:np){

      for(f in 1:nf){
        # calculate pstar quantile of TAC recommendation dist
        TACused[,p,f] <- apply(Data_p_A[[p]][[f]]@TAC, 2, quantile, p = MOM@pstar, na.rm = T) #Data_p_A[[p]][[f]]@TAC#TAC_A[,p,f]#apply(Data_p_A[[p]][[f]]@TAC, 2, quantile, p = MOM@pstar, na.rm = T)

        checkNA[p,f,y] <- sum(is.na(TACused[,p,f]))

        # LastEi[,p,f] <- rep(1,nsim) # no effort adjustment
        LastTAE[,p,f] <-  rep(NA, nsim) # no current TAE exists
        histTAE[,p,f] <- rep(NA, nsim) # no existing TAE
        LastSpatial[,p,f,] <- array(MPA[p,f,nyears,], dim=c(nareas, nsim)) #
        LastAllocat[,p,f] <- rep(1, nsim) # default assumption of reallocation of effort to open areas
        LastCatch[,p,f] <- apply(CB[,p,f,,nyears,], 1, sum)

        Effort_pot[,p,f] <- rep(NA, nsim) # No bio-economic model

        MPCalcs <-  DLMtool::CalcMPDynamics(MPRecs=MPRecs_A[[p]][[f]], y=y, nyears=nyears, proyears=proyears, nsim=nsim,
                                           LastTAE=LastTAE[,p,f], histTAE=histTAE[,p,f],
                                           LastSpatial=LastSpatial[,p,f,], LastAllocat=LastAllocat[,p,f], LastTAC=LastCatch[,p,f],
                                  TACused=TACused[,p,f], maxF=maxF,
                                  LR5_P=FleetPars[[p]][[f]]$LR5_P, LFR_P=FleetPars[[p]][[f]]$LFR_P, Rmaxlen_P=FleetPars[[p]][[f]]$Rmaxlen_P,
                                  retL_P=FleetPars[[p]][[f]]$retL_P, retA_P=FleetPars[[p]][[f]]$retA_P,
                                  L5_P=FleetPars[[p]][[f]]$L5_P, LFS_P=FleetPars[[p]][[f]]$LFS_P, Vmaxlen_P=FleetPars[[p]][[f]]$Vmaxlen_P,
                                  SLarray_P=FleetPars[[p]][[f]]$SLarray_P, V_P=FleetPars[[p]][[f]]$V_P,
                                  Fdisc_P=StockPars[[p]]$Fdisc_P, DR_P=FleetPars[[p]][[f]]$DR_P,
                                  M_ageArray=StockPars[[p]]$M_ageArray,
                                  FM_P=FleetPars[[p]][[f]]$FM_P, FM_Pret=FleetPars[[p]][[f]]$FM_Pret,
                                  Z_P=FleetPars[[p]][[f]]$Z_P, CB_P=FleetPars[[p]][[f]]$CB_P, CB_Pret=FleetPars[[p]][[f]]$CB_Pret,
                                  TAC_f=ImpPars[[p]][[f]]$TAC_f, E_f=ImpPars[[p]][[f]]$E_f, SizeLim_f=ImpPars[[p]][[f]]$SizeLim_f,
                                  VBiomass_P=StockPars[[p]]$VBiomass_P, Biomass_P=StockPars[[p]]$Biomass_P, FinF=FleetPars[[p]][[f]]$FinF,
                                  Spat_targ=FleetPars[[p]][[f]]$Spat_targ,
                                  CAL_binsmid=StockPars[[p]]$CAL_binsmid, Linf=StockPars[[p]]$Linf, Len_age=StockPars[[p]]$Len_age,
                                  maxage=StockPars[[p]]$maxage, nareas=StockPars[[p]]$nareas, Asize=StockPars[[p]]$Asize,
                                  nCALbins=StockPars[[p]]$nCALbins,
                                  qs=FleetPars[[p]][[f]]$qs, qvar=FleetPars[[p]][[f]]$qvar, qinc=FleetPars[[p]][[f]]$qinc,
                                  Effort_pot=Effort_pot[,p,f])

        if(length(SexPars)>0) MPCalcs<-MPCalcsNAs(MPCalcs) # Zeros caused by SexPars

        TACa[,p,f, mm, y] <- TACused[,p,f]#MPCalcs$TACrec # recommended TAC
        LastSpatial[,p,f,] <- MPCalcs$Si
        LastAllocat[,p,f] <- MPCalcs$Ai

        LastTAE[,p,f] <- MPCalcs$TAE # TAE set by MP
        LastCatch[,p,f] <- MPCalcs$TACrec # TAC et by MP

        Effort[,p,f, mm, y] <- rep(MPCalcs$Effort,nsim)[1:nsim]
        FleetPars[[p]][[f]]$CB_P <- MPCalcs$CB_P # removals
        FleetPars[[p]][[f]]$CB_Pret <- MPCalcs$CB_Pret # retained catch
        FleetPars[[p]][[f]]$FM_P <- MPCalcs$FM_P # fishing mortality
        FM_P[,p,f,,,]<- MPCalcs$FM_P

        FleetPars[[p]][[f]]$FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality
        FMret_P[,p,f,,,]<- MPCalcs$FM_Pret
        #FretA[,p,f,,]<- MPCalcs$FM_Pret
        FleetPars[[p]][[f]]$Z_P <- MPCalcs$Z_P # total mortality
        FleetPars[[p]][[f]]$retA_P <- MPCalcs$retA_P # retained-at-age

        FleetPars[[p]][[f]]$retL_P <- MPCalcs$retL_P # retained-at-length
        FleetPars[[p]][[f]]$V_P <- MPCalcs$V_P  # vulnerable-at-age
        VF[,p,f,,]<- MPCalcs$V_P
        FleetPars[[p]][[f]]$SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length

      }
    }

    upyrs <- 1 + (0:(floor(proyears/interval[mm]) - 1)) * interval[mm]  # the years in which there are updates (every three years)
    #upyrs <- 1 + (0:(floor(proyears/2) - 1)) * 2
    if(!silent) {
      cat(".")
      flush.console()
    }

    # --------------------------------------------------------------------------------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------------------------------------------------------------------------------

    # --- Begin projection years ----
    for (y in 2:proyears) {

      # y<-y+1
      if(!silent) {
        cat(".")
        flush.console()
      }

      # -- Calculate MSY stats for this year ----
      if (AnnualMSY) { #
        for(p in 1:np){

          for(f in 1:nf){
            V_Pt[,f,,]<-FleetPars[[p]][[f]]$V_P*apply(CB[,p,f,,nyears,], 1, sum) # Weighted by catch frac
          }
          V_P<-nlz(apply(V_Pt,c(1,3,4),sum),c(1,3),"max") #summed over fleets and normalized to 1

          MSYrefsYr <- sapply(1:nsim, DLMtool::optMSY_eq, StockPars[[p]]$M_ageArray,  StockPars[[p]]$Wt_age,  StockPars[[p]]$Mat_age,
                                V_P, StockPars[[p]]$maxage, StockPars[[p]]$R0, StockPars[[p]]$SRrel, StockPars[[p]]$hs, yr.ind=(nyears+y)+(-1:0))
          MSY_P[,p,mm,y] <- MSYrefsYr[1, ]
          FMSY_P[,p,mm,y] <- MSYrefsYr[2,]
          SSBMSY_P[,p,mm,y] <- MSYrefsYr[3,]
        }

      } # end of annual MSY

      TACa[,,, mm, y] <- TACa[,,, mm, y-1] # TAC same as last year unless changed

      SAYRt <- as.matrix(expand.grid(1:nsim, 1:maxage, y + nyears, 1:nareas))  # Trajectory year
      SAYt <- SAYRt[, 1:3]
      SAYtMP <- cbind(SAYt, mm)
      SYt <- SAYRt[, c(1, 3)]
      SAY1R <- as.matrix(expand.grid(1:nsim, 1:maxage, y - 1, 1:nareas))
      SAYR <- as.matrix(expand.grid(1:nsim, 1:maxage, y, 1:nareas))
      SY <- SAYR[, c(1, 3)]
      SA <- SAYR[, 1:2]
      S1 <- SAYR[, 1]

      SAY <- SAYR[, 1:3]
      S <- SAYR[, 1]
      SR <- SAYR[, c(1, 4)]
      SA2YR <- as.matrix(expand.grid(1:nsim, 2:maxage, y, 1:nareas))
      SA1YR <- as.matrix(expand.grid(1:nsim, 1:(maxage - 1), y -1, 1:nareas))

      for(p in 1:np){
        Perr[,p]<-StockPars[[p]]$Perr_y[,y+nyears+maxage-1]
        M_agecur_y[,p,]<-StockPars[[p]]$M_ageArray[,,nyears+y]
        Mat_agecur_y[,p,]<-StockPars[[p]]$Mat_age[,,nyears+y]
        K[,p]<-StockPars[[p]]$Karray[,nyears+y]
        Linf[,p]<-StockPars[[p]]$Linfarray[,nyears+y]
        t0[,p]<-StockPars[[p]]$t0array[,nyears+y]
      }

      NextYrN <- sapply(1:nsim, function(x)
           popdynOneMICE(np,nf,nareas, maxage, Ncur=array(N_P[x,,,y-1,],c(np,maxage,nareas)), Vcur=array(VF[x,,,,nyears+y-1],c(np,nf,maxage)),
                    #Retcur=array(FretA[x,,,,nyears+y-1],c(np,nf,maxage)),
                    #Fcur=apply(apply(array(FM_P[x,,,,y-1,],c(np,nf,maxage,nareas)),c(1,2,4),max),c(1,2),mean),   # note that Fcur is apical F but, in popdynOneMICE it is DIVIDED in future years between the two areas depending on vulnerabile biomass. So to get Fcur you need to sum over areas (a bit weird)
                    FMretx=array(FMret_P[x,,,,y-1,],c(np,nf,maxage,nareas)),
                    FMx=array(FM_P[x,,,,y-1,],c(np,nf,maxage,nareas)),   # note that Fcur is apical F but, in popdynOneMICE it is DIVIDED in future years between the two areas depending on vulnerabile biomass. So to get Fcur you need to sum over areas (a bit weird)
                    PerrYrp=Perr[x,], hsx=hs[x,], aRx=matrix(aR[x,,],nrow=np), bRx=matrix(bR[x,,],nrow=np),
                    movy=array(mov[x,,,,,nyears+y],c(np,maxage,nareas,nareas)), Spat_targ=array(Spat_targ_y[x,,],c(np,nf)), SRrelx=SRrel[x,],
                    M_agecur=matrix(M_agecur_y[x,,],nrow=np), Mat_agecur=matrix(Mat_agecur_y[x,,],nrow=np),
                    Asizex=matrix(Asize[x,,],ncol=nareas), Kx =K[x,], Linfx=Linf[x,],t0x=t0[x,],Mx=M[x,],
                    R0x=R0[x,],R0ax=matrix(R0a[x,,],nrow=np),SSBpRx=matrix(SSBpR[x,,],nrow=np),ax=a_y,
                    bx=b_y,Rel=Rel,SexPars=SexPars,x=x))

      # x=1; Ncur=array(N_P[x,,,y-1,],c(np,maxage,nareas)); Vcur=array(VF[x,,,,nyears+y-1],c(np,nf,maxage)); FMretx=array(FMret_P[x,,,,y-1,],c(np,nf,maxage,nareas));  FMx=array(FM_P[x,,,,y-1,],c(np,nf,maxage,nareas));   PerrYrp=Perr[x,]; hsx=hs[x,]; aRx=matrix(aR[x,,],nrow=np); bRx=matrix(bR[x,,],nrow=np);  movy=array(mov[x,,,,,nyears+y],c(np,maxage,nareas,nareas)); Spat_targ=array(Spat_targ_y[x,,],c(np,nf)); SRrelx=SRrel[x,]; M_agecur=matrix(M_agecur_y[x,,],nrow=np); Mat_agecur=matrix(Mat_agecur_y[x,,],nrow=np);   Asizex=matrix(Asize[x,,],ncol=nareas); Kx =K[x,]; Linfx=Linf[x,];t0x=t0[x,];Mx=M[x,]; R0x=R0[x,]; R0ax=matrix(R0a[x,,],nrow=np); SSBpRx=matrix(SSBpR[x,,],nrow=np); ax=a_y; bx=b_y; Rel=Rel; SexPars=SexPars; x=x

      N_P[,,,y,]<-aperm(array(as.numeric(unlist(NextYrN[1,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
      Biomass_P[,,,y,]<-aperm(array(as.numeric(unlist(NextYrN[23,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
      SSN_P[,,,y,]<-aperm(array(as.numeric(unlist(NextYrN[24,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
      SSB_P[,,,y,]<-aperm(array(as.numeric(unlist(NextYrN[25,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))

      VBiomass_P[,,,y,]<-aperm(array(as.numeric(unlist(NextYrN[19,], use.names=FALSE)), dim=c(np,maxage, nareas, nsim)), c(4,1,2,3))
      FML <- apply(array(FM_P[, ,,, y-1, ],c(nsim,np,nf,maxage,nareas)), c(1, 3), max)

      for(p in 1:np){
        StockPars[[p]]$N_P<-N_P[,p,,,]
        StockPars[[p]]$Biomass_P<-Biomass_P[,p,,,]
        StockPars[[p]]$SSN_P<-SSN_P[,p,,,]
        StockPars[[p]]$SSB_P<-SSB_P[,p,,,]
        StockPars[[p]]$VBiomass_P<-VBiomass_P[,p,,,]
        #for(f in 1:nf)FleetPars[[p]][[f]]$FML<-FML[]
      }

      # --- An update year ----
      if (y %in% upyrs) {
        # rewrite the DLM object and run the TAC function
        yind <- upyrs[match(y, upyrs) - 1]:(upyrs[match(y, upyrs)] - 1)

        for(p in 1:np){

          for(f in 1:nf) V_Pt[,f,,]<-FleetPars[[p]][[f]]$V_P*LastCatch[,p,f] # Weighted by catch frac
          V_P<-nlz(apply(V_Pt,c(1,3,4),sum),c(1,3),"max") #summed over fleets and normalized to 1 # vulnerability
          MLbin <- (StockPars[[p]]$CAL_bins[1:(length(StockPars[[p]]$CAL_bins) - 1)] + StockPars[[p]]$CAL_bins[2:length(StockPars[[p]]$CAL_bins)])/2

          for(f in 1:nf){

            # use the retained catch
            CBtemp <-  FleetPars[[p]][[f]]$CB_Pret[, , yind, , drop=FALSE] # retained catch-at-age
            CNtemp <- FleetPars[[p]][[f]]$retA_P[,,yind+nyears, drop=FALSE] * apply(StockPars[[p]]$N_P[,,yind,, drop=FALSE], c(1,2,3), sum) # retained age structure

            CBtemp[is.na(CBtemp)] <- tiny
            CBtemp[!is.finite(CBtemp)] <- tiny
            CNtemp[is.na(CNtemp)] <- tiny
            CNtemp[!is.finite(CNtemp)] <- tiny

            Cobs <- Cbiasa[,p,f, nyears + yind] * Cerr[, p,f,nyears + yind] * apply(CBtemp, c(1, 3), sum, na.rm = T)
            Cobs[is.na(Cobs)] <- tiny
            Recobs <- ObsPars[[p]][[f]]$Recerr[, nyears + yind] * apply(array(N_P[,p, 1, yind, ], c(nsim, interval[mm], nareas)), c(1, 2), sum)

            CAA <- array(0, dim = c(nsim, interval[mm], maxage))  # Catch  at age array

            # # a multinomial observation model for catch-at-age data
            for (i in 1:nsim) {
              for (j in 1:interval[mm]) {
                # if (all(CNtemp[i, , j]<1)) { # this is a fix for low sample sizes. If CAA is zero across the board a single fish is caught in age class of model selectivity (dumb I know)
                #   CNtemp[i, floor(maxage/3), j] <- 1
                # }
                CAA[i, j, ] <- ceiling(-0.5 + rmultinom(1, ObsPars[[p]][[f]]$CAA_ESS[i], CNtemp[i, , j]) * ObsPars[[p]][[f]]$CAA_nsamp[i]/ObsPars[[p]][[f]]$CAA_ESS[i])   # a multinomial observation model for catch-at-age data
              }
            }

            ## Calculate CAL ####
            CAL <- array(NA, dim = c(nsim, interval[mm], StockPars[[p]]$nCALbins))  # the catch at length array
            # # a multinomial observation model for catch-at-length data

            vn <- (apply(N_P[,p,,,], c(1,2,3), sum) * FleetPars[[p]][[f]]$retA_P[,,(nyears+1):(nyears+proyears)]) # numbers at age that would be retained
            vn <- aperm(vn, c(1,3,2))

            nyrs <- length(yind)
            tempSize <- lapply(1:nsim, DLMtool::genSizeCompWrap, vn[,yind,, drop=FALSE], CAL_binsmid=StockPars[[p]]$CAL_binsmid, FleetPars[[p]][[f]]$retL_P,
                               ObsPars[[p]][[f]]$CAL_ESS, ObsPars[[p]][[f]]$CAL_nsamp,
                               StockPars[[p]]$Linfarray[,nyears + yind, drop=FALSE],
                               StockPars[[p]]$Karray[,nyears + yind, drop=FALSE],
                               StockPars[[p]]$t0array[,nyears + yind,drop=FALSE], StockPars[[p]]$LenCV, truncSD=2)
            CAL <- aperm(array(as.numeric(unlist(tempSize, use.names=FALSE)), dim=c(length(yind), length(StockPars[[p]]$CAL_binsmid), nsim)), c(3,1,2))

            # calculate LFC - approx 5th percentile
            LFC <- unlist(lapply(tempSize, function(x) DLMtool::getfifth(x[nrow(x), ], StockPars[[p]]$CAL_binsmid)))
            LFC[is.na(LFC)] <- 1
            LFC[LFC<1] <- 1

            I2 <- (cbind(apply(Biomass[,p,,,], c(1, 3), sum), apply(Biomass_P[,p,,,], c(1, 3), sum)[, 1:(y - 1)])^ObsPars[[p]][[f]]$betas) * ObsPars[[p]][[f]]$Ierr[, 1:(nyears + (y - 1))]
            I2[is.na(I2)] <- tiny
            I2 <- I2/apply(I2, 1, mean)

            Depletion <- apply(StockPars[[p]]$SSB_P[, , y, ], 1, sum)/StockPars[[p]]$SSB0 # apply(SSB[, , 1, ], 1, sum)
            Depletion[Depletion < tiny] <- tiny

            NextYrNtemp <- lapply(1:nsim, function(x)
              popdynOneTScpp(nareas, maxage, SSBcurr=colSums(SSB_P[x,p,,y, ]), Ncurr=N_P[x,p,,y,],
                             Zcurr=matrix(StockPars[[p]]$M_ageArray[x,,y+nyears], ncol=nareas, nrow=maxage),
                             PerrYr=StockPars[[p]]$Perr_y[x, y+nyears+maxage-1], hs=StockPars[[p]]$hs[x],
                             R0a=StockPars[[p]]$R0a[x,], SSBpR=StockPars[[p]]$SSBpR[x,], aR=StockPars[[p]]$aR[x,], bR=StockPars[[p]]$bR[x,],
                             mov=StockPars[[p]]$mov[x,,,,nyears+y], SRrel=StockPars[[p]]$SRrel[x]))

            N_PNext <- aperm(array(unlist(NextYrNtemp), dim=c(maxage, nareas, nsim, 1)), c(3,1,4,2))
            VBiomassNext <- VBiomass_P[,p,,,]
            SPAYt<-cbind(SAYt[,1],rep(p,nrow(SAYt)),SAYt[,2:3])
            VBiomassNext[SAYR] <- N_PNext * StockPars[[p]]$Wt_age[SAYt] * V_P[SAYt]  # Calculate vulnerable for abundance

            A <- apply(VBiomassNext[, , y, ], 1, sum)
            # A <- apply(VBiomass_P[, , y, ], 1, sum)

            A[is.na(A)] <- tiny
            Asp <- apply(SSB_P[,p, , y, ], 1, sum)  # SSB Abundance
            Asp[is.na(Asp)] <- tiny
            OFLreal <- A * FMSY_P[,p,mm,y]

            # - update data object ----
            # assign all the new data
            MSElist[[p]][[f]][[mm]]@OM$A <- A
            MSElist[[p]][[f]][[mm]]@Year <- 1:(nyears + y - 1)
            MSElist[[p]][[f]][[mm]]@Cat <- cbind(MSElist[[p]][[f]][[mm]]@Cat, Cobs)
            MSElist[[p]][[f]][[mm]]@Ind <- I2
            MSElist[[p]][[f]][[mm]]@Rec <- cbind(MSElist[[p]][[f]][[mm]]@Rec, Recobs)
            MSElist[[p]][[f]][[mm]]@t <- rep(nyears + y, nsim)
            MSElist[[p]][[f]][[mm]]@AvC <- apply(MSElist[[p]][[f]][[mm]]@Cat, 1, mean)
            MSElist[[p]][[f]][[mm]]@Dt <- ObsPars[[p]][[f]]$Dbias * Depletion * rlnorm(nsim, mconv(1, ObsPars[[p]][[f]]$Derr), sdconv(1, ObsPars[[p]][[f]]$Derr))
            oldCAA <- MSElist[[p]][[f]][[mm]]@CAA
            MSElist[[p]][[f]][[mm]]@CAA <- array(0, dim = c(nsim, nyears + y - 1, maxage))
            MSElist[[p]][[f]][[mm]]@CAA[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAA[, 1:(nyears + y - interval[mm] - 1), ] # there is some bug here sometimes oldCAA (MSElist[[p]][[f]][[mm]]@CAA previously) has too many years of observations
            MSElist[[p]][[f]][[mm]]@CAA[, nyears + yind, ] <- CAA
            MSElist[[p]][[f]][[mm]]@Dep <- ObsPars[[p]][[f]]$Dbias * Depletion * rlnorm(nsim, mconv(1, ObsPars[[p]][[f]]$Derr), sdconv(1, ObsPars[[p]][[f]]$Derr))
            MSElist[[p]][[f]][[mm]]@Abun <- A * ObsPars[[p]][[f]]$Abias * rlnorm(nsim, mconv(1, ObsPars[[p]][[f]]$Aerr), sdconv(1, ObsPars[[p]][[f]]$Aerr))
            MSElist[[p]][[f]][[mm]]@SpAbun <- Asp * ObsPars[[p]][[f]]$Abias * rlnorm(nsim, mconv(1, ObsPars[[p]][[f]]$Aerr), sdconv(1, ObsPars[[p]][[f]]$Aerr))
            MSElist[[p]][[f]][[mm]]@CAL_bins <- StockPars[[p]]$CAL_bins
            oldCAL <- MSElist[[p]][[f]][[mm]]@CAL
            MSElist[[p]][[f]][[mm]]@CAL <- array(0, dim = c(nsim, nyears + y - 1, StockPars[[p]]$nCALbins))
            MSElist[[p]][[f]][[mm]]@CAL[, 1:(nyears + y - interval[mm] - 1), ] <- oldCAL[, 1:(nyears + y - interval[mm] - 1), ]# there is some bug here: sometimes oldCAL (MSElist[[p]][[f]][[mm]]@CAL previously) has too many years of observations
            MSElist[[p]][[f]][[mm]]@CAL[, nyears + yind, ] <- CAL[, 1:interval[mm], ]

            temp <- CAL * rep(MLbin, each = nsim * interval[mm])
            MSElist[[p]][[f]][[mm]]@ML <- cbind(MSElist[[p]][[f]][[mm]]@ML, apply(temp, 1:2, sum)/apply(CAL, 1:2, sum))
            MSElist[[p]][[f]][[mm]]@Lc <- cbind(MSElist[[p]][[f]][[mm]]@Lc, array(MLbin[apply(CAL, 1:2, which.max)], dim = c(nsim, interval[mm])))
            nuCAL <- CAL
            for (i in 1:nsim) for (j in 1:interval[mm]) nuCAL[i, j, 1:match(max(1, MSElist[[p]][[f]][[mm]]@Lc[i, j]), MLbin,nomatch=1)] <- NA
            temp <- nuCAL * rep(MLbin, each = nsim * interval[mm])
            MSElist[[p]][[f]][[mm]]@Lbar <- cbind(MSElist[[p]][[f]][[mm]]@Lbar, apply(temp,1:2, sum, na.rm=TRUE)/apply(nuCAL, 1:2, sum, na.rm=TRUE))

            MSElist[[p]][[f]][[mm]]@LFC <- LFC * ObsPars[[p]][[f]]$LFCbias
            MSElist[[p]][[f]][[mm]]@LFS <- FleetPars[[p]][[f]]$LFS[nyears + y,] * ObsPars[[p]][[f]]$LFSbias

            # update growth, maturity estimates for current year
            MSElist[[p]][[f]][[mm]]@vbK <-  StockPars[[p]]$Karray[, nyears+y] * ObsPars[[p]][[f]]$Kbias
            MSElist[[p]][[f]][[mm]]@vbt0 <- StockPars[[p]]$t0 * ObsPars[[p]][[f]]$t0bias

            MSElist[[p]][[f]][[mm]]@vbLinf <- StockPars[[p]]$Linfarray[, nyears+y] * ObsPars[[p]][[f]]$Linfbias
            MSElist[[p]][[f]][[mm]]@L50 <- StockPars[[p]]$L50array[, nyears+y] * ObsPars[[p]][[f]]$lenMbias
            MSElist[[p]][[f]][[mm]]@L95 <- StockPars[[p]]$L95array[, nyears+y] * ObsPars[[p]][[f]]$lenMbias
            MSElist[[p]][[f]][[mm]]@L95[is.na(MSElist[[p]][[f]][[mm]]@L95)]<-MSElist[[p]][[f]][[mm]]@vbLinf # this is just to robustify 'numbers models' like Grey Seal that do not generate (and will never use) real length observations
            MSElist[[p]][[f]][[mm]]@L95[MSElist[[p]][[f]][[mm]]@L95 > 0.9 * MSElist[[p]][[f]][[mm]]@vbLinf] <- 0.9 * MSElist[[p]][[f]][[mm]]@vbLinf[MSElist[[p]][[f]][[mm]]@L95 > 0.9 * MSElist[[p]][[f]][[mm]]@vbLinf]  # Set a hard limit on ratio of L95 to Linf
            MSElist[[p]][[f]][[mm]]@L50[MSElist[[p]][[f]][[mm]]@L50 > 0.9 * MSElist[[p]][[f]][[mm]]@L95] <- 0.9 * MSElist[[p]][[f]][[mm]]@L95[MSElist[[p]][[f]][[mm]]@L50 > 0.9 * MSElist[[p]][[f]][[mm]]@L95]  # Set a hard limit on ratio of L95 to Linf

            MSElist[[p]][[f]][[mm]]@Ref <- OFLreal
            MSElist[[p]][[f]][[mm]]@Ref_type <- "Simulated OFL"
            MSElist[[p]][[f]][[mm]]@Misc <- Data_p_A[[p]][[f]]@Misc
            MSElist[[p]][[f]][[mm]]@MPrec <- TACa[,p,f, mm, y] # last MP  TAC recommendation
            MSElist[[p]][[f]][[mm]]@MPeff <- Effort[,p,f, mm, y-1] # last recommended effort

          } # end of fleet
        } # end of stock

        # assign('Data',DataList[[mm]],envir=.GlobalEnv) # for debugging fun

        if(MPcond=="MMP"){
          # returns a hierarchical list object stock then fleet of MPrec objects

          DataList<-getDataList(MSElist,mm) # returns a hierarchical list object stock then fleet of Data objects
          MPRecs_A <- applyMMP(DataList, MP = MPs[mm], reps = 1, silent=TRUE)  # # returns a hierarchical list object stock then fleet then slot type of Rec
          Data_p_A <- MPrecs_A_blank
          for(p in 1:np)for(f in 1:nf){
            Data_p_A[[p]][[f]]<-MSElist[[p]][[f]][[mm]]
            Data_p_A[[p]][[f]]@TAC<-MPRecs_A[[p]][[f]]$TAC # record TAC rec in Data
          }

        }else if(MPcond=="complex"){

          MPRecs_A <- Data_p_A <- MPrecs_A_blank # A temporary blank hierarchical list object stock by fleet
          realVB<-abind::abind(apply(VBiomass[,,,1:nyears,],c(1,2,4),sum,na.rm=T),apply(VBiomass_P[,,,1:(y-1),],c(1,2,4),sum,na.rm=T),along=3) # need this for aggregating data and distributing TACs over stocks

          curdat<-multiDataS(MSElist,StockPars,np,mm,nf,realVB)
          runMP <- applyMP(curdat, MPs = MPs[mm], reps = 1, silent=TRUE)  # Apply MP

          Stock_Alloc<-realVB[,,nyears]/apply(realVB[,,nyears],1,sum)

          for(p in 1:np)for(f in 1:nf){

            MPRecs_A[[p]][[f]]<-runMP[[1]][[1]]
            MPRecs_A[[p]][[f]]$TAC<-runMP[[1]][[1]]$TAC*MOM@Allocation[[p]][,f]*Stock_Alloc[,p]
            MPRecs_A[[p]][[f]]$Effort<-runMP[[1]][[1]]$Effort*MOM@Efactor[[p]][,f]

            if(length(MPRecs_A[[p]][[f]]$Effort)>0)if(is.na(MPRecs_A[[p]][[f]]$Effort[1,1])) MPRecs_A[[p]][[f]]$Effort<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$Effort))
            if(length(MPRecs_A[[p]][[f]]$TAC)>0)if(is.na(MPRecs_A[[p]][[f]]$TAC[1,1])) MPRecs_A[[p]][[f]]$TAC<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))
            if(is.na(MPRecs_A[[p]][[f]]$Spatial[1,1])) MPRecs_A[[p]][[f]]$Spatial<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))

            Data_p_A[[p]][[f]]<-runMP[[2]]
            Data_p_A[[p]][[f]]@TAC<-MPRecs_A[[p]][[f]]$TAC

          }

        }else{

          MPRecs_A <- Data_p_A <- MPrecs_A_blank # A temporary blank hierarchical list object stock by fleet

          for(p in 1:np){

            if(MPcond=="bystock"){

              if(nf>1){
                curdat<-multiData(MSElist,StockPars,p,mm,nf)
              }else{
                curdat<-MSElist[[p]][[f]][[mm]]
              }

              runMP <- applyMP(curdat, MPs = MPs[[p]][mm], reps = MOM@reps, silent=TRUE)  # Apply MP

              # Do allocation calcs
              TAC_A[,p,]<-array(as.vector(unlist(runMP[[1]][[1]]$TAC))*MOM@Allocation[[p]],c(nsim,nf))
              TAE_A[,p,]<-array(as.vector(unlist(runMP[[1]][[1]]$Effort))*MOM@Efactor[[p]],c(nsim,nf))

              for(f in 1:nf){
                MPRecs_A[[p]][[f]]<-runMP[[1]][[1]]
                MPRecs_A[[p]][[f]]$TAC<-matrix(TAC_A[,p,f],nrow=1) # Just pass the allocated TAC
                MPRecs_A[[p]][[f]]$Effort<-matrix(TAE_A[,p,f],nrow=1)
                # This next line is to make the NULL effort recommendations of an output control MP compatible with CalcMPdynamics (expects a null matrix)
                if(is.na(MPRecs_A[[p]][[f]]$Effort[1,1])) MPRecs_A[[p]][[f]]$Effort<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$Effort))
                if(is.na(MPRecs_A[[p]][[f]]$TAC[1,1])) MPRecs_A[[p]][[f]]$TAC<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))
                if(is.na(MPRecs_A[[p]][[f]]$Spatial[1,1])) MPRecs_A[[p]][[f]]$Spatial<-matrix(NA,nrow=0,ncol=ncol(MPRecs_A[[p]][[f]]$TAC))
                Data_p_A[[p]][[f]]<-runMP[[2]]
                Data_p_A[[p]][[f]]@TAC<-MPRecs_A[[p]][[f]]$TAC # Copy the allocated TAC
              }

            }else if(MPcond=="byfleet"){

              for(f in 1:nf){

                curdat<-MSElist[[p]][[f]][[mm]]
                runMP <- applyMP(curdat, MPs = MPrefs[p,f,mm], reps = MOM@reps, silent=TRUE)  # Apply MP
                MPRecs_A[[p]][[f]]<-runMP[[1]][[1]]
                Data_p_A[[p]][[f]]<-runMP[[2]]
                Data_p_A[[p]][[f]]@TAC <- MPRecs_A[[p]][[f]]$TAC

              }

            }

          } # end of stocks

        } # end of MMP?

        for(p in 1:np){

          for(f in 1:nf){

            # calculate pstar quantile of TAC recommendation dist
            TACused[,p,f] <- apply(Data_p_A[[p]][[f]]@TAC, 2, quantile, p = MOM@pstar, na.rm = T)

            checkNA[p,f,y] <-checkNA[p,f,y] + sum(is.na(TACused[,p,f]))


            MPCalcs <- DLMtool::CalcMPDynamics(MPRecs=MPRecs_A[[p]][[f]], y=y, nyears=nyears, proyears=proyears, nsim=nsim,
                                               LastTAE=LastTAE[,p,f], histTAE=histTAE[,p,f],
                                               LastSpatial=LastSpatial[,p,f,], LastAllocat=LastAllocat[,p,f], LastTAC=LastCatch[,p,f],
                                               TACused=TACused[,p,f], maxF=maxF,
                                               LR5_P=FleetPars[[p]][[f]]$LR5_P, LFR_P=FleetPars[[p]][[f]]$LFR_P, Rmaxlen_P=FleetPars[[p]][[f]]$Rmaxlen_P,
                                               retL_P=FleetPars[[p]][[f]]$retL_P, retA_P=FleetPars[[p]][[f]]$retA_P,
                                               L5_P=FleetPars[[p]][[f]]$L5_P, LFS_P=FleetPars[[p]][[f]]$LFS_P, Vmaxlen_P=FleetPars[[p]][[f]]$Vmaxlen_P,
                                               SLarray_P=FleetPars[[p]][[f]]$SLarray_P, V_P=FleetPars[[p]][[f]]$V_P,
                                               Fdisc_P=StockPars[[p]]$Fdisc_P, DR_P=FleetPars[[p]][[f]]$DR_P,
                                               M_ageArray=StockPars[[p]]$M_ageArray,
                                               FM_P=FleetPars[[p]][[f]]$FM_P, FM_Pret=FleetPars[[p]][[f]]$FM_Pret,
                                               Z_P=FleetPars[[p]][[f]]$Z_P, CB_P=FleetPars[[p]][[f]]$CB_P, CB_Pret=FleetPars[[p]][[f]]$CB_Pret,
                                               TAC_f=ImpPars[[p]][[f]]$TAC_f, E_f=ImpPars[[p]][[f]]$E_f, SizeLim_f=ImpPars[[p]][[f]]$SizeLim_f,
                                               VBiomass_P=StockPars[[p]]$VBiomass_P, Biomass_P=StockPars[[p]]$Biomass_P, FinF=FleetPars[[p]][[f]]$FinF,
                                               Spat_targ=FleetPars[[p]][[f]]$Spat_targ,
                                               CAL_binsmid=StockPars[[p]]$CAL_binsmid, Linf=StockPars[[p]]$Linf, Len_age=StockPars[[p]]$Len_age,
                                               maxage=StockPars[[p]]$maxage, nareas=StockPars[[p]]$nareas, Asize=StockPars[[p]]$Asize,
                                               nCALbins=StockPars[[p]]$nCALbins,
                                               qs=FleetPars[[p]][[f]]$qs, qvar=FleetPars[[p]][[f]]$qvar, qinc=FleetPars[[p]][[f]]$qinc,
                                               Effort_pot=Effort_pot[,p,f])

            if(length(SexPars)>0) MPCalcs<-MPCalcsNAs(MPCalcs) # Zeros caused by SexPars

            TACa[,p,f, mm, y] <- MPCalcs$TACrec # recommended TAC
            LastSpatial[,p,f,] <- MPCalcs$Si
            LastAllocat[,p,f] <- MPCalcs$Ai

            LastTAE[,p,f] <- MPCalcs$TAE # adjustment to TAE
            LastCatch[,p,f] <- MPCalcs$TACrec

            Effort[,p,f, mm, y] <- rep(MPCalcs$Effort,nsim)[1:nsim]
            FleetPars[[p]][[f]]$CB_P <- MPCalcs$CB_P # removals
            FleetPars[[p]][[f]]$CB_Pret <- MPCalcs$CB_Pret # retained catch
            FleetPars[[p]][[f]]$FM_P <- MPCalcs$FM_P # fishing mortality
            FM_P[,p,f,,,]<- MPCalcs$FM_P
            FleetPars[[p]][[f]]$FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality
            # FretA[,p,f,,]<- MPCalcs$FM_Pret
            FleetPars[[p]][[f]]$Z_P <- MPCalcs$Z_P # total mortality
            FleetPars[[p]][[f]]$retA_P <- MPCalcs$retA_P # retained-at-age

            FleetPars[[p]][[f]]$retL_P <- MPCalcs$retL_P # retained-at-length
            FleetPars[[p]][[f]]$V_P <- MPCalcs$V_P  # vulnerable-at-age
            VF[,p,f,,]<- MPCalcs$V_P
            FleetPars[[p]][[f]]$SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length

            # apply combined MP ----
            if("DataOut"%in%names(control))if(control$DataOut == y) return(MSElist)

            # calculate pstar quantile of TAC recommendation dist
          } # end of fleets
        } # end of stocks


      } else {
        # --- Not an update yr ----

        for(p in 1:np){

          for(f in 1:nf){

            NoMPRecs <- MPRecs_A[[p]][[f]]
            # NoMPRecs[lapply(NoMPRecs, length) > 0 ] <- NULL
            NoMPRecs$Spatial <- NA

            MPCalcs <- DLMtool::CalcMPDynamics(MPRecs=NoMPRecs, y=y, nyears=nyears, proyears=proyears, nsim=nsim,
                                               LastTAE=LastTAE[,p,f], histTAE=histTAE[,p,f],
                                               LastSpatial=LastSpatial[,p,f,], LastAllocat=LastAllocat[,p,f], LastTAC=LastCatch[,p,f],
                                      TACused=TACused[,p,f], maxF=maxF,
                                      LR5_P=FleetPars[[p]][[f]]$LR5_P, LFR_P=FleetPars[[p]][[f]]$LFR_P, Rmaxlen_P=FleetPars[[p]][[f]]$Rmaxlen_P,
                                      retL_P=FleetPars[[p]][[f]]$retL_P, retA_P=FleetPars[[p]][[f]]$retA_P,
                                      L5_P=FleetPars[[p]][[f]]$L5_P, LFS_P=FleetPars[[p]][[f]]$LFS_P, Vmaxlen_P=FleetPars[[p]][[f]]$Vmaxlen_P,
                                      SLarray_P=FleetPars[[p]][[f]]$SLarray_P, V_P=FleetPars[[p]][[f]]$V_P,
                                      Fdisc_P=StockPars[[p]]$Fdisc_P, DR_P=FleetPars[[p]][[f]]$DR_P,
                                      M_ageArray=StockPars[[p]]$M_ageArray,
                                      FM_P=FleetPars[[p]][[f]]$FM_P, FM_Pret=FleetPars[[p]][[f]]$FM_Pret,
                                      Z_P=FleetPars[[p]][[f]]$Z_P, CB_P=FleetPars[[p]][[f]]$CB_P, CB_Pret=FleetPars[[p]][[f]]$CB_Pret,
                                      TAC_f=ImpPars[[p]][[f]]$TAC_f, E_f=ImpPars[[p]][[f]]$E_f, SizeLim_f=ImpPars[[p]][[f]]$SizeLim_f,
                                      VBiomass_P=StockPars[[p]]$VBiomass_P, Biomass_P=StockPars[[p]]$Biomass_P, FinF=FleetPars[[p]][[f]]$FinF,
                                      Spat_targ=FleetPars[[p]][[f]]$Spat_targ,
                                      CAL_binsmid=StockPars[[p]]$CAL_binsmid, Linf=StockPars[[p]]$Linf, Len_age=StockPars[[p]]$Len_age,
                                      maxage=StockPars[[p]]$maxage, nareas=StockPars[[p]]$nareas, Asize=StockPars[[p]]$Asize,
                                      nCALbins=StockPars[[p]]$nCALbins,
                                      qs=FleetPars[[p]][[f]]$qs, qvar=FleetPars[[p]][[f]]$qvar, qinc=FleetPars[[p]][[f]]$qinc,
                                      Effort_pot=Effort_pot[,p,f])

            if(length(SexPars)>0) MPCalcs<-MPCalcsNAs(MPCalcs) # Zeros caused by SexPars

            TACa[,p,f, mm, y] <- TACused[,p,f] # recommended TAC
            #TACa[,p,f, mm, y] <- MPCalcs$TACrec # recommended TAC
            LastSpatial[,p,f,] <- MPCalcs$Si
            LastAllocat[,p,f] <- MPCalcs$Ai

            LastTAE[,p,f] <- MPCalcs$TAE
            # LastEi[,p,f] <- MPCalcs$Ei # adjustment to effort
            LastCatch[,p,f] <- MPCalcs$TACrec

            Effort[,p,f, mm, y] <- rep(MPCalcs$Effort,nsim)[1:nsim]
            FleetPars[[p]][[f]]$CB_P <- MPCalcs$CB_P # removals
            FleetPars[[p]][[f]]$CB_Pret <- MPCalcs$CB_Pret # retained catch
            FleetPars[[p]][[f]]$FM_P <- MPCalcs$FM_P # fishing mortality
            FM_P[,p,f,,,]<- MPCalcs$FM_P
            FleetPars[[p]][[f]]$FM_Pret <- MPCalcs$FM_Pret # retained fishing mortality
            FMret_P[,p,f,,,]<- MPCalcs$FM_Pret
            #FretA[,p,f,,]<- MPCalcs$FM_Pret
            FleetPars[[p]][[f]]$Z_P <- MPCalcs$Z_P # total mortality
            FleetPars[[p]][[f]]$retA_P <- MPCalcs$retA_P # retained-at-age

            FleetPars[[p]][[f]]$retL_P <- MPCalcs$retL_P # retained-at-length
            FleetPars[[p]][[f]]$V_P <- MPCalcs$V_P  # vulnerable-at-age
            VF[,p,f,,]<- MPCalcs$V_P
            FleetPars[[p]][[f]]$SLarray_P <- MPCalcs$SLarray_P # vulnerable-at-length

          }  # end of fleets

        } # end of stocks

      }  # end of not an update year

      #checkNA[y] <- sum(is.na(TACused))

    }  # end of year

    B_BMSYa[, ,mm, ] <- apply(SSB_P, c(1,2, 4), sum, na.rm=TRUE)/array(SSBMSY_P[,,mm,],c(nsim,np,proyears))  # SSB relative to SSBMSY

    for(p in 1:np)for(f in 1:nf)FMa[,p,f, mm, ] <- -log(1 - apply(FleetPars[[p]][[f]]$CB_P, c(1, 3), sum, na.rm=TRUE)/apply(VBiomass_P[,p,,,]+FleetPars[[p]][[f]]$CB_P, c(1, 3), sum, na.rm=TRUE))
    for(f in 1:nf)F_FMSYa[, ,f,mm, ] <- FMa[,,f, mm, ]/FMSY_P[,,mm,]

    Ba[, ,mm, ] <- apply(Biomass_P, c(1, 2,4), sum, na.rm=TRUE) # biomass
    SSBa[, ,mm, ] <- apply(SSB_P, c(1, 2,4), sum, na.rm=TRUE) # spawning stock biomass
    VBa[, ,mm, ] <- apply(VBiomass_P, c(1, 2, 4), sum, na.rm=TRUE) # vulnerable biomass

    for(p in 1:np)for(f in 1:nf)Ca[, p,f,mm, ] <- apply(FleetPars[[p]][[f]]$CB_P, c(1, 3), sum, na.rm=TRUE) # removed
    for(p in 1:np)for(f in 1:nf)CaRet[, p,f,mm, ] <- apply(FleetPars[[p]][[f]]$CB_Pret, c(1, 3), sum, na.rm=TRUE) # retained catch

    # Store Pop and Catch-at-age and at-length for last projection year
    PAAout[ ,, mm, ] <- apply(array(N_P[ , , , proyears, ],c(nsim,np,maxage,nareas)), c(1,2,3), sum) # population-at-age

    for(p in 1:np){for(f in 1:nf){
      CNtemp <- apply(FleetPars[[p]][[f]]$CB_Pret, c(1,2,3), sum)/StockPars[[p]]$Wt_age[(nyears+1):nyears+proyears]
      CAAout[ ,p,f, mm, ] <- CNtemp[,,proyears] # nsim, maxage # catch-at-age
      lastyr<-dim(MSElist[[p]][[f]][[mm]]@CAL)[2]
      CALout[[p]][[f]][[mm]] <- MSElist[[p]][[f]][[mm]]@CAL[,lastyr,] # catch-at-length in last year
    }}

    if (!silent) {
      cat("\n")
      if (all(checkNA != nsim) & !all(checkNA == 0)) {
        # print number of NAs
        # message(checkNA)
        # message(checkNA[upyrs])
        ntot <- sum(checkNA[,,upyrs])
        totyrs <- sum(checkNA[,,upyrs] >0)
        nfrac <- round(ntot/(length(upyrs)*nsim),2)*100

        message(totyrs, ' years had TAC = NA for some simulations (', nfrac, "% of total simulations)")
        message('Used TAC_y = TAC_y-1')
      }

    }

    if (!parallel)
      if("progress"%in%names(control))
        if(control$progress)
          shiny::incProgress(1/nMP, detail = round(mm*100/nMP))

  }  # end of mm methods

  # Miscellaneous reporting
  if(PPD)Misc<-MSElist

  # rescale effort to today
  for(p in 1:np)for(f in 1:nf)Effort[,p,f,,] <-  array(Effort[,p,f,,],c(nsim,nMP,proyears))/array(FleetPars[[p]][[f]]$FinF, dim=c(nsim, nMP, proyears))
  OM<-Obsout<-CALbins<-list()
  for(p in 1:np){
    OM[[p]]<-Obsout[[p]]<-list()
    CALbins[[p]]<-StockPars[[p]]$CAL_binsmid
    for(f in 1:nf){
      OM[[p]][[f]]<-DataList[[p]][[f]]@OM
      Obsout[[p]][[f]]<-DataList[[p]][[f]]@Obs
    }
  }

  Misc[['MOM']]<-MOM

  if(class(MPs)=="character")MPs<-list(MPs) # need to reformat MMP and complex mode to work with MSEout slot

  ## Create MSE Object ####
  MSEout <- new("MMSE", Name = MOM@Name, nyears, proyears, nMPs=nMP, MPs=MPs, MPcond=MPcond,MPrefs=MPrefs,nsim, nstocks=np, nfleets=nf,
                Snames=Snames, Fnames=Fnames, Stocks=Stocks, Fleets=Fleets, Obss=Obs, Imps=Imps,OM=OM, Obs=Obsout, B_BMSY=B_BMSYa, F_FMSY=F_FMSYa, B=Ba,
                SSB=SSBa, VB=VBa, FM=FMa, CaRet, TAC=TACa, SSB_hist = SSB, CB_hist = CB,
                FM_hist = FM, Effort = Effort, PAA=PAAout, CAA=CAAout, CAL=CALout, CALbins=CALbins,
                MSY_P = MSY_P, FMSY_P = FMSY_P, SSBMSY_P = SSBMSY_P,   Misc = Misc)


  # Store MSE info
  attr(MSEout, "version") <- paste0("MSEtool: ",packageVersion("MSEtool"), "  DLMtool: ",packageVersion("DLMtool"))
  attr(MSEout, "date") <- date()
  attr(MSEout, "R.version") <- R.version

  MSEout

}




