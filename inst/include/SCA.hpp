
#ifndef SCA_hpp
#define SCA_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type SCA(objective_function<Type> *obj) {

  using namespace ns_SCA;

  DATA_VECTOR(C_hist);    // Total catch
  DATA_SCALAR(rescale);   // Scalar for R0
  DATA_VECTOR(I_hist);    // Index
  DATA_MATRIX(CAA_hist);  // Catch-at-age proportions
  DATA_VECTOR(CAA_n);     // Annual samples in CAA
  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_VECTOR(M);         // Natural mortality at age
  DATA_VECTOR(weight);    // Weight-at-age at the beginning of the year
  DATA_VECTOR(mat);       // Maturity-at-age at the beginning of the year
  DATA_STRING(vul_type);  // String indicating whether logistic or dome vul is used
  DATA_STRING(I_type);    // String whether index surveys B, VB, or SSB
  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_STRING(CAA_dist);  // String indicating whether CAA is multinomial or lognormal
  DATA_IVECTOR(est_early_rec_dev);
  DATA_IVECTOR(est_rec_dev); // Indicator of whether rec_dev is estimated in model or fixed at zero
  DATA_INTEGER(yindF);    // Year for which to estimate F, all other F are deviations from this F

  PARAMETER(R0x);
  PARAMETER(transformed_h);
  PARAMETER(F_equilibrium);
  PARAMETER_VECTOR(vul_par);

  PARAMETER_VECTOR(logF);

  PARAMETER(log_sigma);
  PARAMETER(log_omega);
  PARAMETER(log_tau);
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

  Type sigma = exp(log_sigma);
  Type omega = exp(log_omega);
  Type tau = exp(log_tau);

  Type penalty = 0;
  Type prior = 0.;

  // Vulnerability
  vector<Type> vul(max_age);
  if(vul_type == "logistic") {
    vul = calc_logistic_vul(vul_par, max_age, prior);
  } else {
    vul = calc_dome_vul(vul_par, max_age, prior);
  }

  ////// Equilibrium reference points and per-recruit quantities
  vector<Type> NPR_virgin = calc_NPR(Type(0), vul, M, max_age);

  Type EPR0 = sum_EPR(NPR_virgin, weight, mat);

  Type B0 = R0 * sum_BPR(NPR_virgin, weight);
  Type N0 = R0 * NPR_virgin.sum();
  Type E0 = R0 * EPR0;
  Type VB0 = R0 * sum_VBPR(NPR_virgin, weight, vul);

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
  Type Arec = CR/EPR0;
  Brec /= E0;

  ////// During time series year = 1, 2, ..., n_y
  matrix<Type> N(n_y+1, max_age);   // Numbers at year and age
  matrix<Type> CAApred(n_y, max_age);   // Catch (in numbers) at year and age at the mid-point of the season
  vector<Type> CN(n_y);             // Catch in numbers
  vector<Type> Cpred(n_y);
  vector<Type> F(n_y);
  vector<Type> Ipred(n_y);          // Predicted index at year
  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(max_age-1);
  vector<Type> VB(n_y+1);           // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year

  CN.setZero();
  Cpred.setZero();
  VB.setZero();
  B.setZero();
  E.setZero();

  // Equilibrium quantities (leading into first year of model)
  vector<Type> NPR_equilibrium = calc_NPR(F_equilibrium, vul, M, max_age);
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
  for(int a=0;a<max_age;a++) {
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
  F(yindF) = exp(logF(yindF));
  for(int y=0;y<n_y;y++) {
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y), h, R0, E0);
    } else {
      R(y+1) = Ricker_SR(E(y), h, R0, E0);
    }

    if(y<n_y-1) {
      if(est_rec_dev(y+1)) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * tau * tau);
    }
    N(y+1,0) = R(y+1);

    if(y != yindF) {
      Type Ftmp = F(yindF) * exp(logF(y));
      Type tmp2 = 3 - Ftmp;
      F(y) = CppAD::CondExpLt(tmp2, Type(0), 3 - posfun(tmp2, Type(0), penalty), Ftmp);
    }

    for(int a=0;a<max_age;a++) {
      Type meanN = N(y,a) * (1 - exp(-vul(a) * F(y) - M(a))) / (vul(a) * F(y) + M(a));
      CAApred(y,a) = vul(a) * F(y) * meanN;
      if(a<max_age-1) N(y+1,a+1) = N(y,a) * exp(-vul(a) * F(y) - M(a));
      if(a==max_age-1) N(y+1,a) += N(y,a) * exp(-vul(a) * F(y) - M(a));
      CN(y) += CAApred(y,a);
      Cpred(y) += CAApred(y,a) * weight(a);
      B(y+1) += N(y+1,a) * weight(a);
      VB(y+1) += N(y+1,a) * weight(a) * vul(a);
      E(y+1) += N(y+1,a) * weight(a) * mat(a);
    }
  }

  // Calculate nuisance parameters and likelihood
  Type q;
  if(I_type == "B") {
    q = calc_q(I_hist, B);
    for(int y=0;y<n_y;y++) Ipred(y) = q * B(y);
  } else if(I_type == "VB") {
    q = calc_q(I_hist, VB);
    for(int y=0;y<n_y;y++) Ipred(y) = q * VB(y);
  } else {
    q = calc_q(I_hist, E);
    for(int y=0;y<n_y;y++) Ipred(y) = q * E(y);
  }

  vector<Type> nll_comp(4);
  nll_comp.setZero();
  for(int y=0;y<n_y;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll_comp(0) -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
    if(C_hist(y) > 0) {
      if(!R_IsNA(asDouble(CAA_n(y)))) {
        vector<Type> loglike_CAAobs(max_age);
        vector<Type> loglike_CAApred(max_age);
        loglike_CAApred = CAApred.row(y)/CN(y);
        loglike_CAAobs = CAA_hist.row(y);
        if(CAA_dist == "multinomial") {
          loglike_CAAobs *= CAA_n(y);
          nll_comp(1) -= dmultinom_(loglike_CAAobs, loglike_CAApred, true);
        } else {
          nll_comp(1) -= dlnorm_comp(loglike_CAAobs, loglike_CAApred);
        }
      }
      nll_comp(2) -= dnorm(log(C_hist(y)), log(Cpred(y)), omega, true);
    }
    if(est_rec_dev(y)) nll_comp(3) -= dnorm(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<max_age-1;a++) {
    if(est_early_rec_dev(a)) nll_comp(3) -= dnorm(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_comp.sum() + penalty + prior;

  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(omega);
  ADREPORT(sigma);
  ADREPORT(tau);
  ADREPORT(q);

  REPORT(omega);
  REPORT(sigma);
  REPORT(tau);

  REPORT(NPR_virgin);
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
