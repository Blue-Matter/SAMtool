
#ifndef cDD_hpp
#define cDD_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type cDD(objective_function<Type> *obj) {
  using namespace ns_cDD;
  
  DATA_SCALAR(Winf);
  DATA_SCALAR(Kappa);
  DATA_INTEGER(ny);
  DATA_INTEGER(k);
  DATA_SCALAR(wk);
  DATA_VECTOR(C_hist);
  DATA_SCALAR(dep);
  DATA_SCALAR(rescale);
  DATA_MATRIX(I_hist);
  DATA_IVECTOR(I_units);
  DATA_MATRIX(I_sd);
  DATA_VECTOR(MW_hist);
  DATA_STRING(SR_type);
  DATA_INTEGER(n_itF);
  DATA_VECTOR(LWT);
  DATA_INTEGER(nsurvey);
  DATA_INTEGER(fix_sigma);
  DATA_INTEGER(state_space);
  DATA_IVECTOR(use_prior); // Boolean vector, whether to set a prior for R0, h, M, q (length of 3 + nsurvey)
  DATA_MATRIX(prior_dist); // Distribution of priors for R0, h, M, q (rows), columns indicate parameters of distribution calculated in R (see RCM_prior fn)

  PARAMETER(R0x);
  PARAMETER(transformed_h);
  PARAMETER(log_M);
  PARAMETER(F_equilibrium);
  PARAMETER(log_sigma);
  PARAMETER(log_sigma_W);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_rec_dev);

  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else h = exp(transformed_h);
  h += 0.2;
  Type R0 = exp(R0x)/rescale;
  Type M = exp(log_M);
  Type sigma = exp(log_sigma);
  Type sigma_W = exp(log_sigma_W);
  Type tau = exp(log_tau);
  int SR_type2 = SR_type == "BH";


  //--DECLARING DERIVED VALUES
  Type BPR0 = cDD_BPR(Type(0), M, wk, Kappa, Winf);
  Type B0 = BPR0 * R0;
  Type N0 = R0/M;

  Type CR, Brec;

  if(SR_type == "BH") {
    CR = 4 *h;
    CR /= 1-h;
    Brec = 5*h - 1;
    Brec /= (1-h) * B0;
  } else {
    CR = pow(5*h, 1.25);
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= B0;
  }
  Type Arec = CR/BPR0;

  //--DECLARING STORAGE VECTORS
  vector<Type> B(ny+1);
  vector<Type> N(ny+1);
  vector<Type> R(ny+k);

  vector<Type> F(ny);
  vector<Type> Z(ny);
  vector<Type> Cpred(ny);
  matrix<Type> Ipred(ny,nsurvey);
  vector<Type> MWpred(ny);

  vector<Type> BPRinf(ny);
  vector<Type> Rec_dev(ny);
  vector<Type> Binf(ny);
  vector<Type> Ninf(ny);

  //--INITIALIZE
  Type BPReq = cDD_BPR(F_equilibrium, M, wk, Kappa, Winf);
  Type Req = cDD_R(BPReq, Arec, Brec, SR_type2);

  B(0) = Req * BPReq;
  N(0) = Req/(F_equilibrium + M);
  for(int tt=0;tt<k;tt++) R(tt) = Req;
  if(state_space) {
    Rec_dev(0) = exp(log_rec_dev(0) - 0.5 * tau * tau);
    SIMULATE { 
      Rec_dev(0) = exp(rnorm(log_rec_dev(0) - 0.5 * tau * tau, tau));
    }
    R(0) *= Rec_dev(0);
  } else {
    Rec_dev.fill(1);
  }
  
  Type Ceqpred = F_equilibrium * B(0);

  Type penalty = 0; // Penalty to likelihood for high F > max_F
  Type prior = 0;
  if(dep > 0) { // Penalty for initial depletion to get corresponding F
    prior = -dnorm_(log(B(0)/B0), log(dep), Type(0.01), true);
  }
  
  for(int tt=0; tt<ny; tt++) {
    Type F_start = CppAD::CondExpLe(C_hist(tt), Type(1e-8), Type(0), -log(1 - C_hist(tt)/B(tt)));
    F(tt) = cDD_F(F_start, C_hist(tt), M, Winf, Kappa, wk, N, B, Cpred, BPRinf, Binf, R, Ninf,
      CppAD::Integer(CppAD::CondExpLe(C_hist(tt), Type(1e-8), Type(1), Type(n_itF))), tt);
    Z(tt) = F(tt) + M;
    MWpred(tt) = B(tt)/N(tt);

    N(tt+1) = Ninf(tt) + (N(tt) - Ninf(tt)) * exp(-Z(tt));
    B(tt+1) = Binf(tt);
    B(tt+1) += Kappa * Winf * (N(tt) - Ninf(tt)) / (Z(tt) + Kappa);
    B(tt+1) += (B(tt) - Binf(tt) - Kappa * Winf * (N(tt) - Ninf(tt))/(Z(tt) + Kappa)) * exp(-Z(tt) - Kappa);

    if(SR_type == "BH") {
      R(tt+k) = BH_SR(B(tt), h, R0, B0);
    } else {
      R(tt+k) = Ricker_SR(B(tt), h, R0, B0);
    }
    if(state_space && tt<ny-1) {
      Rec_dev(tt+1) = exp(log_rec_dev(tt+1) - 0.5 * tau * tau);
      SIMULATE {
        Rec_dev(tt+1) = exp(rnorm(log_rec_dev(tt+1) - 0.5 * tau * tau, tau));
      }
      R(tt+1) *= Rec_dev(tt+1);
    }
    
    SIMULATE {
      C_hist(tt) = Cpred(tt);
    }
    
  }

  //--ARGUMENTS FOR NLL
  vector<Type> q = calc_q(I_hist, B, N, Ipred, nsurvey, I_units, ny);

  // Objective function
  //creates storage for jnll and sets value to 0
  vector<Type> nll_comp(nsurvey+2);
  nll_comp.setZero();

  for(int sur=0;sur<nsurvey;sur++) {
    for(int tt=0;tt<ny;tt++) {
      if(LWT(sur) > 0 && !R_IsNA(asDouble(I_hist(tt,sur)))) {
        if(fix_sigma) {
          nll_comp(sur) -= dnorm_(log(I_hist(tt,sur)), log(Ipred(tt,sur)), I_sd(tt,sur), true);
          
          SIMULATE {
            I_hist(tt,sur) = exp(rnorm(log(Ipred(tt,sur)), I_sd(tt,sur)));
          }
          
        } else {
          nll_comp(sur) -= dnorm(log(I_hist(tt,sur)), log(Ipred(tt,sur)), sigma, true);
          
          SIMULATE {
            I_hist(tt,sur) = exp(rnorm(log(Ipred(tt,sur)), sigma));
          }
        }
      }
    }
    nll_comp(sur) *= LWT(sur);
  }
  
  for(int tt=0;tt<ny;tt++) {
    if(LWT(nsurvey) > 0 && !R_IsNA(asDouble(MW_hist(tt)))) {
      nll_comp(nsurvey) -= dnorm(log(MW_hist(tt)), log(MWpred(tt)), sigma_W, true);
      
      SIMULATE {
        MW_hist(tt) = exp(rnorm(log(MWpred(tt)), sigma_W));
      }
    }
  }
  nll_comp(nsurvey) *= LWT(nsurvey);
  if(state_space) {
    for(int tt=0;tt<ny;tt++) nll_comp(nsurvey+1) -= dnorm(log_rec_dev(tt), Type(0), tau, true);
  }

  //Summing individual jnll and penalties
  prior -= calc_prior(use_prior, prior_dist, R0x, h, SR_type == "BH", log_M, q, rescale);
  Type nll = nll_comp.sum() + penalty + prior;

  //-------REPORTING-------//
  ADREPORT(R0);
  ADREPORT(h);
  if(CppAD::Variable(log_M)) ADREPORT(M);
  ADREPORT(q);
  if(CppAD::Variable(log_sigma)) ADREPORT(sigma);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(MW_hist.sum() > 0) {
    ADREPORT(sigma_W);
    REPORT(sigma_W);
  }
  if(!fix_sigma) REPORT(sigma);
  if(state_space) REPORT(tau);
  REPORT(nll);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(CR);
  REPORT(BPR0);
  REPORT(Cpred);
  REPORT(Ceqpred);
  REPORT(Ipred);
  REPORT(MWpred);
  REPORT(q);
  REPORT(B);
  REPORT(N);
  REPORT(R);
  REPORT(F);
  REPORT(M);
  REPORT(Z);
  REPORT(BPRinf);
  REPORT(Binf);
  REPORT(Ninf);
  REPORT(h);
  REPORT(R0);
  REPORT(N0);
  REPORT(B0);
  REPORT(log_rec_dev);
  REPORT(Rec_dev);
  REPORT(nll_comp);
  REPORT(penalty);
  REPORT(prior);
  
  SIMULATE {
    REPORT(C_hist);
    REPORT(I_hist);
    REPORT(MW_hist);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
