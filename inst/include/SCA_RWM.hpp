
#ifndef SCA_RWM_hpp
#define SCA_RWM_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type SCA_RWM(objective_function<Type> *obj) {

  using namespace ns_SCA;

  DATA_VECTOR(C_hist);    // Total catch
  DATA_SCALAR(rescale);   // Scalar for R0
  DATA_MATRIX(I_hist);    // Index
  DATA_MATRIX(I_sd);
  DATA_IVECTOR(I_units);
  DATA_MATRIX(I_vul);
  DATA_IVECTOR(abs_I);
  DATA_INTEGER(nsurvey);
  DATA_VECTOR(LWT);
  DATA_MATRIX(CAA_hist);  // Catch-at-age proportions
  DATA_VECTOR(CAA_n);     // Annual samples in CAA
  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(n_age);    // Maximum age (plus-group)
  //DATA_VECTOR(M);       // Natural mortality at age
  DATA_VECTOR(weight);    // Weight-at-age at the beginning of the year
  DATA_VECTOR(mat);       // Maturity-at-age at the beginning of the year
  DATA_STRING(vul_type);  // String indicating whether logistic or dome vul is used
  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_STRING(CAA_dist);  // String indicating whether CAA is multinomial or lognormal
  DATA_IVECTOR(est_early_rec_dev);
  DATA_IVECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero
  DATA_INTEGER(yindF);    // Year for which to estimate F, all other F are deviations from this F
  DATA_VECTOR(M_bounds);  // Lower and upper bounds for M walk

  PARAMETER(R0x);
  PARAMETER(transformed_h);
  PARAMETER(F_equilibrium);
  PARAMETER_VECTOR(vul_par);

  PARAMETER_VECTOR(log_F_dev);
  
  PARAMETER(log_M_start);
  PARAMETER_VECTOR(logit_M_walk);

  PARAMETER(log_omega);
  PARAMETER(log_tau);
  PARAMETER(log_tau_M);
  
  PARAMETER_VECTOR(log_early_rec_dev);
  PARAMETER_VECTOR(log_rec_dev);

  Type R0 = exp(R0x)/rescale;
  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else {
    h = exp(transformed_h);
  }
  h += 0.2;

  Type omega = exp(log_omega);
  Type tau = exp(log_tau);
  Type tau_M = exp(log_tau_M);
  
  Type penalty = 0;
  Type prior = 0.;

  // Vulnerability
  vector<Type> vul(n_age);
  if(vul_type == "logistic") {
    vul = calc_logistic_vul(vul_par, n_age, prior);
  } else {
    vul = calc_dome_vul(vul_par, n_age, prior);
  }
  
  // Set up M vector
  vector<Type> M(n_y);
  vector<Type> logit_M(n_y);
  vector< vector<Type>> NPR0(n_y);
  vector<Type> EPR0(n_y);
  
  M(0) = exp(log_M_start);
  logit_M(0) = logit2(M(0), M_bounds(0), M_bounds(1), M(0));
  for(int y=0;y<n_y;y++) {
    if(y > 0) {
      logit_M(y) = logit_M(y-1) + logit_M_walk(y-1);
      M(y) = invlogit2(logit_M(y), M_bounds(0), M_bounds(1), M(0));
    }
    NPR0(y) = calc_NPR(Type(0), vul, M(y), n_age);
    EPR0(y) = sum_EPR(NPR0(y), weight, mat);
  }
  
  ////// Equilibrium reference points and per-recruit quantities
  Type B0 = R0 * sum_BPR(NPR0(0), weight);
  Type N0 = R0 * NPR0(0).sum();
  Type E0 = R0 * EPR0(0);
  Type VB0 = R0 * sum_VBPR(NPR0(0), weight, vul);

  Type CR, Brec;
  if(SR_type == "BH") {
    CR = 4 * h;
    CR /= 1-h;
    Brec = 5*h - 1;
    Brec /= (1-h);
  } else {
    CR = pow(5*h, 1.25);
    Brec = 1.25 * log(5*h);
  }
  Type Arec = CR/EPR0(0);
  Brec /= E0;

  ////// During time series year = 1, 2, ..., n_y
  matrix<Type> N(n_y+1, n_age);   // Numbers at year and age
  matrix<Type> CAApred(n_y, n_age);   // Catch (in numbers) at year and age at the mid-point of the season
  vector<Type> CN(n_y);             // Catch in numbers
  vector<Type> Cpred(n_y);
  vector<Type> F(n_y);
  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(n_age-1);
  vector<Type> VB(n_y+1);           // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year
  
  CN.setZero();
  Cpred.setZero();
  VB.setZero();
  B.setZero();
  E.setZero();
  
  // Equilibrium quantities (leading into first year of model)
  vector<Type> NPR_equilibrium = calc_NPR(F_equilibrium, vul, M(0), n_age);
  Type EPR_eq = sum_EPR(NPR_equilibrium, weight, mat);
  Type R_eq;

  if(SR_type == "BH") {
    R_eq = Arec * EPR_eq - 1;
  } else {
    R_eq = log(Arec * EPR_eq);
  }
  R_eq /= Brec * EPR_eq;

  R(0) = R_eq;
  if(est_rec_dev(0)) R(0) *= exp(log_rec_dev(0) - 0.5 * tau * tau);
  for(int a=0;a<n_age;a++) {
    if(a == 0) {
      N(0,a) = R(0) * NPR_equilibrium(a);
    } else {
      R_early(a-1) = R_eq;
      if(est_early_rec_dev(a-1)) R_early(a-1) *= exp(log_early_rec_dev(a-1) - 0.5 * tau * tau);
      N(0,a) = R_early(a-1) * NPR_equilibrium(a);
    }
    B(0) += N(0,a) * weight(a);
    VB(0) += N(0,a) * weight(a) * vul(a);
    E(0) += N(0,a) * weight(a) * mat(a);
  }

  // Loop over all other years
  F(yindF) = exp(log_F_dev(yindF));
  for(int y=0;y<n_y;y++) {
    // Calculate this year's F
    if(y != yindF) {
      Type Ftmp = F(yindF) * exp(log_F_dev(y));
      F(y) = CppAD::CondExpLt(3 - Ftmp, Type(0), 3 - posfun(3 - Ftmp, Type(0), penalty), Ftmp);
    }
    
    // Calculate this year's catch, CAA, and next year's abundance and SSB (ex. age-0)
    for(int a=0;a<n_age;a++) {
      Type meanN = N(y,a) * (1 - exp(-vul(a) * F(y) - M(y))) / (vul(a) * F(y) + M(y));
      CAApred(y,a) = vul(a) * F(y) * meanN;
      if(a<n_age-1) {
        N(y+1,a+1) = N(y,a) * exp(-vul(a) * F(y) - M(y));
      } else {
        N(y+1,a) += N(y,a) * exp(-vul(a) * F(y) - M(y));
      }
      CN(y) += CAApred(y,a);
      Cpred(y) += CAApred(y,a) * weight(a);
      E(y+1) += N(y+1,a) * weight(a) * mat(a);
    }
    
    // Calculate next year's recruitment, total biomass, and vulnerable biomass
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y+1), h, R0, E0);
    } else {
      R(y+1) = Ricker_SR(E(y+1), h, R0, E0);
    }
    if(y<n_y-1) {
      if(est_rec_dev(y+1)) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * tau * tau);
    }
    N(y+1,0) = R(y+1);

    for(int a=0;a<n_age;a++) {
      B(y+1) += N(y+1,a) * weight(a);
      VB(y+1) += N(y+1,a) * weight(a) * vul(a);
    }
  }

  // Calculate nuisance parameters and likelihood
  // Ipred updated in calc_q function
  vector<Type> q(nsurvey);
  array<Type> s_CAA(n_y,n_age,nsurvey);
  matrix<Type> s_CN(n_y,nsurvey);
  matrix<Type> s_BN(n_y,nsurvey);
  matrix<Type> Ipred(n_y,nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      for(int a=0;a<n_age;a++) {
        if(I_vul.col(sur).sum() > 0) {
          s_CAA(y,a,sur) = I_vul(a,sur) * N(y,a);
        } else {
          s_CAA(y,a,sur) = vul(y,a) * N(y,a);
        }
        s_CN(y,sur) += s_CAA(y,a,sur);
        if(I_units(sur)) s_BN(y,sur) += s_CAA(y,a,sur) * weight(a); // Biomass vulnerable to survey
      }
    }
    if(!I_units(sur)) s_BN.col(sur) = s_CN.col(sur); // Abundance vulnerable to survey
    q(sur) = calc_q(I_hist, s_BN, sur, sur, Ipred, abs_I, n_y); // This function updates Ipred
  }
  
  vector<Type> nll_comp(nsurvey+3);
  nll_comp.setZero();
  for(int y=0;y<n_y;y++) {
    for(int sur=0;sur<nsurvey;sur++) {
      if(LWT(sur) > 0 && !R_IsNA(asDouble(I_hist(y,sur)))) {
        nll_comp(sur) -= LWT(sur) * dnorm(log(I_hist(y,sur)), log(Ipred(y,sur)), I_sd(y,sur), true);
      }
    }
    if(C_hist(y) > 0) {
      if(!R_IsNA(asDouble(CAA_n(y)))) {
        vector<Type> loglike_CAAobs(n_age);
        vector<Type> loglike_CAApred(n_age);
        loglike_CAApred = CAApred.row(y)/CN(y);
        loglike_CAAobs = CAA_hist.row(y);
        if(CAA_dist == "multinomial") {
          loglike_CAAobs *= CAA_n(y);
          nll_comp(nsurvey) -= dmultinom_(loglike_CAAobs, loglike_CAApred, true);
        } else {
          nll_comp(nsurvey) -= dlnorm_comp(loglike_CAAobs, loglike_CAApred);
        }
      }
      nll_comp(nsurvey+1) -= dnorm(log(C_hist(y)), log(Cpred(y)), omega, true);
    }
    if(est_rec_dev(y)) nll_comp(nsurvey+2) -= dnorm(log_rec_dev(y), Type(0), tau, true);
    if(y<n_y-1) nll_comp(nsurvey+2) -= dnorm(logit_M_walk(y), Type(0), tau_M, true);
  }
  for(int a=0;a<n_age-1;a++) {
    if(est_early_rec_dev(a)) nll_comp(nsurvey+2) -= dnorm(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_comp.sum() + penalty + prior;

  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(omega);
  ADREPORT(tau);
  ADREPORT(tau_M);
  ADREPORT(q);
  ADREPORT(logit_M);

  REPORT(omega);
  REPORT(tau);
  REPORT(tau_M);
  
  REPORT(logit_M);
  REPORT(logit_M_walk);
  REPORT(M);
  
  REPORT(NPR0);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(EPR0);
  REPORT(CR);
  REPORT(h);
  REPORT(R0);
  REPORT(B0);
  REPORT(N0);
  REPORT(E0);
  REPORT(VB0);

  REPORT(vul_par);
  REPORT(vul);

  REPORT(F);
  REPORT(q);

  REPORT(N);
  REPORT(CN);
  REPORT(Cpred);
  REPORT(CAApred);
  REPORT(Ipred);
  REPORT(R);
  REPORT(R_early);
  REPORT(VB);
  REPORT(B);
  REPORT(E);

  REPORT(log_early_rec_dev);
  REPORT(log_rec_dev);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
