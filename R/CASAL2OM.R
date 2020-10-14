
# ====================================================================
# === CASAL2OM =======================================================
# ====================================================================

# April 2021

#' Reads MLE estimates from CASAL file structure into an operating model
#'
#'
#' @description A (prototype) function that uses the file location of a fitted CASAL assessment model including input files to population the
#' various slots of an operating model with MLE parameter estimates. The function mainly populates the Stock and Fleet portions
#' of the operating model; the user still needs to parameterize most of the observation and implementation portions of the operating model.
#' @param CASALdir A folder with CASAL input and output files in it
#' @param Obs The observation model (class Obs).
#' @param Imp The implementation model (class Imp).
#' @param Name The common name of the operating model
#' @param Agency The fishery management agency
#' @param Region The geographical location
#' @param Sponsor Who funded the work
#' @param Latitude In degrees north
#' @param Longitude In degrees west
#' @param nsim The number of simulations to take for parameters with uncertainty (for OM@@cpars custom parameters).
#' @param proyears The number of projection years for MSE
#' @param interval The number of years between management updates
#' @param pstar The quantile for TAC management given stochasticity
#' @param maxF The maximum allowable F in the operating model.
#' @param reps The number of stochastic replicates within each simulation in the operating model.
#' @param seed The random seed for the operating model.
#' @param Common_Name The name of the species
#' @param Species The species latin name
#' @param Source Reference to assessment documentation e.g. a url
#' @param Author Who did the assessment
#' @return An object of class OM.
#' @author T. Carruthers
#' @export
#' @seealso \link{SS2OM}
CASAL2OM<-function(CASALdir,Obs=DLMtool::Precise_Unbiased, Imp=DLMtool::Perfect_Imp,
                   Name=NA,Agency=NA, Region=NA, Sponsor=NA, Latitude=NA,
                   Longitude=NA, nsim=48,proyears=50, interval=4, pstar=0.5, maxF=2,
                   reps=1,seed=1, Common_Name=NA,Species=NA, Source=NA,Author=NA){

  # Known issues / areas for improvement
  # - currently MLE only
  # - the maturity ogive ripping only works for type "allvalues_bounded" - user defined minage and maxage assign a vector that long assuming all values below minage are 0 and all above maxage are 1
  # - the mortality rate is only set up here for a single-age invariant M
  # - no sex specific population dynamics are currently included
  # - post release mortality rate was assumed to be tagging mortality

  # For debugging
  # options(max.print=100000);
  # CASALdir="C:/Users/tcarruth/Documents/GitHub/DLMDev/Case_Studies/Z - INCOMPLETE/Toothfish_Patagonia_CCAMLR/data/CASAL/"
  # CASALdir="C:/GitHub/DLMDev/Case_Studies/Z - INCOMPLETE/Toothfish_Patagonia_CCAMLR/data/CASAL/"
  # Obs=DLMtool::Precise_Biased; Imp=DLMtool::Perfect_Imp;  Name=NA;Agency=NA; Region=NA; Sponsor=NA; Latitude=NA;  Longitude=NA; nsim=16;proyears=50; interval=4; pstar=0.5; maxF=2;  reps=1;seed=1;
  # Source="test source";  Common_Name="test CN";Species="test species"; Source="test source";Author="test author"; Name="test name" ;Agency="test agency"; Region="test region"; Sponsor="test sponsor"; Latitude = 50; Longitude=50
  # getlabs<-function(output)sapply(1:length(output),function(x,output)strsplit(output[x]," ")[[1]][1],output=output)
  # getlabs(output)
  # getlabs(mpd)

  # Construct OM
  Stock <- new("Stock")
  Fleet <- new("Fleet")
  OM <- new("OM", Stock = Stock, Fleet = Fleet, Obs = Obs, Imp = Imp)

  # OM Misc
  OM@Name <- Name
  OM@Agency <- Agency
  OM@Region <- Region
  OM@Sponsor <- Sponsor
  OM@Latitude <- Latitude
  OM@Longitude <- Longitude
  OM@nsim <- nsim
  OM@proyears <- proyears
  OM@interval <- interval
  OM@pstar <- pstar
  OM@maxF <- maxF
  OM@reps <- reps
  #cpars
  OM@seed <- seed
  OM@Source <- Source
  OM@Common_Name <- Common_Name
  OM@Species <- Species

  # CASAL pars

  pars<-CASALpars(CASALdir)
  maxage<-pars$maxage
  nyears<-pars$nyears

  OM@maxage<-pars$maxage
  OM@R0<-pars$R0
  OM@M<-rep(pars$M,2)
  OM@Msd<-c(0.005,0.01)
  OM@h<-rep(pars$h,2)
  OM@SRrel<-1 # match(pars$SR_type,c("BH","")) # don't know how CASAL reports the Ricker, maybe just 'R'?
  OM@Perr<-rep(pars$Perr,2)
  #OM@AC
  OM@Linf <-rep(pars$Linf,2)
  OM@K<-rep(pars$K,2)
  OM@t0 <- rep(pars$t0, 2)
  OM@LenCV <-rep(pars$LenCV,2)
  OM@Ksd<-OM@Linfsd<-c(0.01,0.02)

  # Over written by cpars Mat_age
  OM@L50 <-rep(pars$Linf/2,2)
  OM@L50_95 <-rep(pars$Linf/10,2)


  # YOU GOT TO HERE
  Dmu<-pars$SSB[nyears]/pars$SSB[1]
  OM@D<-rep(Dmu,2)
  OM@a<-pars$a
  OM@b<-pars$b
  OM@Size_area_1<-OM@Frac_area_1<-c(0.1,0.1)
  OM@Prob_staying<-c(0.5,0.5)
  OM@Fdisc<-rep(pars$Fdisc,2)
  OM@nyears<-pars$nyears
  OM@Spat_targ<-c(1,1)

  # These are overwritten by cpars Find and V
  OM@EffYears<-c(1,nyears)
  OM@EffLower<-c(0.1,0.99)
  OM@EffUpper<-c(0.11,1)
  OM@Esd <- c(0.1,0.1)
  OM@qinc<-c(0,0)
  OM@qcv<-c(0.01,0.02)
  OM@L5<-rep(pars$Linf/4,2)
  OM@LFS<-rep(pars$Linf/2,2)
  OM@Vmaxlen<-rep(1,2)
  OM@isRel="F"
  #OM@LR5
  #OM@LFR
  #Rmaxlen
  OM@DR<-c(0,0)
  #SelYears
  #AbsSelYears

  recdevs<-log(pars$recdevs)
  procsd<-sd(recdevs)
  procmu <- -0.5 * (procsd)^2  # adjusted log normal mean
  AC<-acf(recdevs,plot=F)$acf[2,1,1]

  OM@Perr<-rep(procsd,2)
  OM@AC=rep(AC,2)

  Perr<-array(NA,c(nsim,nyears+proyears+maxage-1))
  Perr<-matrix(rnorm(nsim*(maxage+nyears+proyears-1),rep(procmu,maxage+nyears+proyears-1),rep(procsd,maxage+nyears+proyears-1)),nrow=nsim)
  Perr[,maxage+(0:(nyears-1))]<-rep(recdevs,each=nsim) # generate a bunch of simulations with uncertainty
  Perr[, 1:(maxage-1)]<-0
  #Perr[]<-0
  #for (y in 2:(maxage-1)) Perr[, y] <- AC * Perr[, y - 1] +   Perr[, y] * (1 - AC * AC)^0.5
  for (y in maxage+(nyears:(nyears + proyears))-1) Perr[, y] <- AC * Perr[, y - 1] +   Perr[, y] * (1 - AC * AC)^0.5
  Perr<-exp(Perr)

  Mat_age <- array(rep(rep(pars$mat_age,each=nsim),pars$nyears+proyears), dim = c(nsim, pars$maxage, pars$nyears + proyears))

  V <- array(NA, dim = c(nsim, maxage, nyears + proyears))
  Vindhist<-as.matrix(expand.grid(1:nsim,1:maxage,1:nyears))
  V[Vindhist]<-pars$VbY[Vindhist[,3:2]]

  # 50 year projection based on average of last three years
  selrecent<-apply(V[1,,nyears+(-2:0)],1,mean) # recent selectivity is mean over last 4 years 2005-2008
  selrecent<-selrecent/max(selrecent)
  Vindpro<-as.matrix(expand.grid(1:nsim,1:maxage,nyears+(1:proyears)))
  V[Vindpro]<-selrecent[Vindpro[,2]]

  Ftrend<-apply(pars$catch,2,sum)/pars$SSB
  Ftrend<-Ftrend/mean(Ftrend)
  Find<-array(rep(Ftrend,each=nsim),c(nsim,nyears))

  OM@CurrentYr<-as.numeric(pars$yearnams[length(pars$yearnams)])

  OM@cpars<-list(V=V,Find=Find,Mat_age=Mat_age,Perr=Perr)

  OM
  #test<-runMSE(OM)

}


#' Rips MLE estimates from CASAL file structure
#'
#'
#' @description A function that uses the file location of a fitted CASAL assessment model including input files
#' to extract data required to populate an OMx class operating model.
#' @param CASALdir A folder with Stock Synthesis input and output files in it
#' @return A list.
#' @author T. Carruthers
#' @export
#' @seealso \link{CASAL2OM}
CASALpars<-function(CASALdir){

  rip<-function(obj,lookup,pos=2,sep=" ",relline=0,cond2=NA){
    if(!is.na(cond2)) obj<-obj[(grep(cond2,obj)+1):length(obj)]
    if(pos!="all"){
      outy<-strsplit(obj[grep(lookup,obj)[1]+relline],sep)[[1]][pos]  # get a value from obj, identified by lookup, in postion pos in the vector, separated by sep, rellines below lookup
    }else{
      outy<-strsplit(obj[grep(lookup,obj)[1]+relline],sep)[[1]]  # get a value from obj, identified by lookup, in postion pos in the vector, separated by sep, rellines below lookup
    }
    outy
  }

  files<-list.files(CASALdir)
  runname<-strsplit(files[grep("casal_estimation",files)],"_")[[1]][1]

  est <- readLines(paste0(CASALdir,"/",files[grep("casal_estimation",files)])) # estimation set up
  out <- readLines(paste0(CASALdir,"/",files[grep("casal_output",files)]))     # what outputs to return
  pop <- readLines(paste0(CASALdir,"/",files[grep("casal_population",files)])) # input data and parameters
  mpd <- readLines(paste0(CASALdir,"/",files[grep("_mpd",files)]))             # q B0 recruitment and selectivity
  output <- readLines(paste0(CASALdir,"/",files[grep("_output.log",files)]))
  #runobj <- readLines(paste0(CASALdir,"/",files[grep("objectives",files)])) # MCMC objective function output
  #runout <- readLines(paste0(CASALdir,"/",files[grep(paste("output",runname,"casal_run",sep="_"),files)])) #
  #runsamp <- readLines(paste0(CASALdir,"/",files[grep("samples",files)]))

  est_type <- rip(est, "@estimator")
  Author   <- rip(output,"User name",sep=": ")
  years    <- rip(out,"years","all")
  years    <- years[2:length(years)]
  nyears   <- length(years)
  yearnams <- rip(out,"years",pos="all")
  yearnams <- yearnams[1+1:nyears]
  maxage   <- as.integer(rip(pop,"@max_age"))
  SR_type  <- rip(pop,"SR")
  Growth_type <- rip(pop,"@size_at_age_type von_Bert")
  R0       <- as.numeric(rip(output,"* R0",pos=1,relline=1))
  Perr     <- as.numeric(rip(pop,"sigma_r"))
  h       <- as.numeric(rip(pop,"steepness"))

  # if all ages
  M        <- as.numeric(rip(pop,"@natural_mortality",relline=1))
  Linf     <- as.numeric(rip(pop,"Linf"))
  t0       <- as.numeric(rip(pop,"t0"))
  K        <- as.numeric(rip(pop,"k"))
  LenCV    <- as.numeric(rip(pop,"cv"))
  a        <- as.numeric(rip(pop,"a "))
  b        <- as.numeric(rip(pop,"b "))
  mat_age<-rep(1,maxage)

  # if @maturity_props
  matinst<- rip(pop,"@maturity_props",pos="all",relline=1)
  mat_age[matinst[3]:matinst[4]]<-matinst[5:length(matinst)]
  mat_age[1:matinst[3]]<-0
  mat_age<-as.numeric(mat_age)

  recdevs<-as.numeric(strsplit(mpd[2]," ")[[1]][3:(2+nyears)])
  SSB<-rip(output,"SSB",pos="all",cond2 = "* SSBs")
  SSB<-as.numeric(SSB[2:length(SSB)])

  Fdisc    <- as.numeric(rip(pop,lookup="@tag ",relline=8))


  # Selectivities
  sel_text<-est[grep("parameter selectivity", est)]
  nf<-length(sel_text)
  selnames<-sapply(1:nf,function(x,sel_text)substr(sel_text[x],start=gregexpr(pattern='\\[',sel_text[x])[[1]][1]+1,stop=gregexpr(pattern='\\]',sel_text[x])[[1]][1]-1),sel_text=sel_text)
  sels<-array(NA,c(nf,maxage))
  for(f in 1:nf)sels[f,]<-as.numeric(rip(output,paste0("selectivity\\[",selnames[f],"\\].all"),pos="all",cond2="Ogive parameter values")[1+1:maxage])

  # CASAL anchors the q to a single fleet and presumably calculates all other qs as nuisance parameters relative to fleet 1 exploitation rate
  fleetnams<-sapply(1:nf,function(x,selnames)strsplit(selnames[x],"_")[[1]][2],selnames=selnames)

  catch<-array(0,c(nf,nyears))
  for(f in 1:nf){
    yrs<-rip(pop,"years", pos="all",cond2=paste0("fishery ",fleetnams[f]))
    yrs<-yrs[2:length(yrs)]
    catches<-rip(pop,"catches", pos="all",cond2=paste0("fishery ",fleetnams[f]))
    catches<-catches[2:length(catches)]
    catch[f,match(yrs,yearnams)]<-as.numeric(catches)
  }

  VbY<-array(NA,c(nyears,maxage))
  for(y in nyears:1){
    Vtemp<-apply(sels*catch[,y],2,sum)
    if(!all(Vtemp==0)){
      VbY[y,]<-Vtemp/max(Vtemp)
    }else{
      VbY[y,]<-VbY[y+1,]
    }
  }

  #for(f in 1:nf)qs[f]<-as.numeric(rip(output,lookup=paste0("q\\[",fleetnams[f],"q\\].q"),pos=1))

  pars <- list(est_type=est_type, Author=Author, years=years,
               nyears=nyears, yearnams=yearnams,maxage=maxage, SR_type=SR_type,Growth_type=Growth_type,
               R0=R0, Perr=Perr, h=h, M=M, Linf=Linf, t0=t0, K=K, LenCV=LenCV,
               a=a, b=b, mat_age=mat_age, Fdisc=Fdisc, recdevs=recdevs, SSB=SSB,
               nf=nf,fleetnams=fleetnams,sels=sels,catch=catch,VbY=VbY)
  pars

}










