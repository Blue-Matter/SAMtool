
#ifndef VPA_hpp
#define VPA_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type VPA(objective_function<Type> *obj) {

  using namespace ns_VPA;

  DATA_VECTOR(I_hist);         // Index
  DATA_MATRIX(CAA_hist);       // Catch-at-age proportions
  DATA_INTEGER(n_y);           // Number of years in model
  DATA_INTEGER(max_age);       // Maximum age (plus-group)
  DATA_VECTOR(M);              // Natural mortality at age
  DATA_VECTOR(weight);         // Weight-at-age at the beginning of the year
  DATA_VECTOR(mat);            // Maturity-at-age at the beginning of the year
  DATA_STRING(vul_type_term);  // Vulnerability function
  DATA_STRING(I_type);         // String whether index surveys B, VB, or SSB
  DATA_INTEGER(nitF);          // The maximum number of iterations to solve for F
  DATA_INTEGER(n_vulpen);      //
  DATA_SCALAR(sigma_vulpen);
  DATA_INTEGER(n_Rpen);
  DATA_SCALAR(sigma_Rpen);

  PARAMETER(logF_term);
  PARAMETER_VECTOR(vul_par);
  PARAMETER(logF_ratio);
  PARAMETER(log_sigma);

  matrix<Type> F(n_y,max_age);
  matrix<Type> N(n_y,max_age);
  matrix<Type> vul(n_y,max_age);
  matrix<Type> CAApred(n_y,max_age);
  vector<Type> Ipred(n_y);

  vector<Type> VB(n_y);
  vector<Type> E(n_y);
  vector<Type> B(n_y);

  E.setZero();
  VB.setZero();
  B.setZero();
  F.setZero();
  N.setZero();

  Type penalty = 0.;
  Type prior = 0.;

  Type F_term = exp(logF_term);
  Type F_ratio = exp(logF_ratio);
  Type sigma = exp(log_sigma);

  // Vulnerability
  // Add option for 'free' vul parameters
  vector<Type> vul_term(max_age);
  if(vul_type_term == "logistic") {
    vul_term = calc_logistic_vul(vul_par, max_age, prior);
  } else if(vul_type_term == "free") {
    for(int a=0;a<max_age-1;a++) vul_term(a) = exp(vul_par(a));
    vul_term(max_age-1) = F_ratio * vul_term(max_age-2);
  } else {
    vul_term = calc_dome_vul(vul_par, max_age, prior);
  }

  // Terminal Year
  for(int a=0;a<max_age;a++) {
    F(n_y-1,a) = vul_term(a) * F_term;
    N(n_y-1,a) = (F(n_y-1,a) + M(a)) * CAA_hist(n_y-1,a);
    N(n_y-1,a) /= (1 - exp(-F(n_y-1,a) - M(a))) * F(n_y-1,a);
  }

  // Backwards recursion of N
  for(int y=n_y-1;y>0;y--) {
    for(int a=1;a<max_age;a++) {
      if(a==max_age-1) {
        F(y-1,a-1) = Newton_VPA_F_plus(F(y,a), F_ratio, M(a-1), M(a), CAA_hist(y-1,a-1), CAA_hist(y-1,a), N(y,a), nitF);
      } else {
        F(y-1,a-1) = CppAD::CondExpGt(CAA_hist(y-1,a-1), Type(1e-4), Newton_VPA_F(F(y,a), M(a-1), CAA_hist(y-1,a-1), N(y,a), nitF), Type(1e-4));
      }

      N(y-1,a-1) = CppAD::CondExpGt(CAA_hist(y-1,a-1), Type(1e-4), (F(y-1,a-1) + M(a-1)) * CAA_hist(y-1,a-1)/(1 - exp(-F(y-1,a-1) - M(a-1)))/F(y-1,a-1), N(y,a) * exp(F(y-1,a-1) + M(a)));

      if(a==max_age-1) {
        F(y-1,a) = F_ratio * F(y-1,a-1);
        N(y-1,a) = (F(y-1,a) + M(a)) * CAA_hist(y-1,a);
        N(y-1,a) /= (1 - exp(-F(y-1,a) - M(a))) * F(y-1,a);
      }
    }
  }

  for(int y=0;y<n_y;y++) {
    vector<Type> Fvec(max_age);
    for(int a=0;a<max_age;a++) {
      CAApred(y,a) = N(y,a) * F(y,a) * (1 - exp(-F(y,a) - M(a)))/(F(y,a) + M(a));
      E(y) += N(y,a) * weight(a) * mat(a);
      B(y) += N(y,a) * weight(a);
      Fvec(a) = F(y,a);
    }
    for(int a=0;a<max_age;a++) {
      vul(y,a) = F(y,a)/max(Fvec);
      VB(y) += N(y,a) * weight(a) * vul(y,a);
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
  //Type sigma = calc_sigma(I_hist, Ipred);

  vector<Type> nll_comp(3);
  nll_comp.setZero();
  for(int y=0;y<n_y;y++) {
    if(!R_IsNA(asDouble(I_hist(y)))) nll_comp(0) -= dnorm(log(I_hist(y)), log(Ipred(y)), sigma, true);
    if(y>n_y-n_vulpen) {
      for(int a=0;a<max_age;a++) {
        nll_comp(1) -= dnorm_(log(vul(y,a)), log(vul(y-1,a)), sigma_vulpen, true);
      }
    }
    if(y>n_y-n_Rpen) nll_comp(2) -= dnorm_(log(N(y,1)), log(N(y-1,1)), sigma_Rpen, true);
  }

  Type fn = nll_comp.sum() + penalty + prior;

  ADREPORT(q);
  ADREPORT(sigma);
  ADREPORT(F_term);
  ADREPORT(F_ratio);

  REPORT(vul_par);
  REPORT(vul_term);
  REPORT(vul);
  REPORT(F);
  REPORT(N);
  REPORT(CAApred);
  REPORT(Ipred);
  REPORT(VB);
  REPORT(E);
  REPORT(B);

  REPORT(q);
  REPORT(sigma);
  REPORT(F_term);
  REPORT(F_ratio);

  REPORT(nll_comp);
  REPORT(penalty);
  REPORT(prior);
  REPORT(fn);

  return fn;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

