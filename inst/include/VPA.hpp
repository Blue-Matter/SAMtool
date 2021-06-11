
#ifndef VPA_hpp
#define VPA_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type VPA(objective_function<Type> *obj) {

  using namespace ns_VPA;

  DATA_MATRIX(I_hist);         // Index
  DATA_MATRIX(I_sd);
  DATA_IVECTOR(I_units);
  DATA_MATRIX(I_vul);
  DATA_IVECTOR(abs_I);
  DATA_INTEGER(nsurvey);
  DATA_VECTOR(LWT);
  DATA_MATRIX(CAA_hist);       // Catch-at-age proportions
  DATA_INTEGER(n_y);           // Number of years in model
  DATA_INTEGER(n_age);       // Maximum age (plus-group)
  DATA_VECTOR(M);              // Natural mortality at age
  DATA_VECTOR(weight);         // Weight-at-age at the beginning of the year
  DATA_STRING(vul_type_term);  // Vulnerability function
  DATA_INTEGER(n_itF);          // The maximum number of iterations to solve for F
  DATA_INTEGER(n_vulpen);
  DATA_SCALAR(vulpen);      //
  DATA_INTEGER(n_Rpen);
  DATA_SCALAR(Rpen);

  PARAMETER(log_Fterm);
  PARAMETER(log_Fratio);
  PARAMETER_VECTOR(vul_par);

  matrix<Type> F_at_age(n_y,n_age);
  matrix<Type> Z_at_age(n_y,n_age);
  vector<Type> F(n_y);
  matrix<Type> N(n_y,n_age);
  matrix<Type> vul(n_y,n_age);
  matrix<Type> CAApred(n_y,n_age);
  vector<Type> VB(n_y);
  vector<Type> B(n_y);
  
  VB.setZero();
  B.setZero();
  F_at_age.setZero();
  Z_at_age.setZero();
  N.setZero();

  Type penalty = 0.;
  Type prior = 0.;

  Type Fterm = exp(log_Fterm);
  Type Fratio = exp(log_Fratio);

  // Vulnerability
  // Add option for 'free' vul parameters
  vector<Type> vul_term(n_age);
  if(vul_type_term == "logistic") {
    vul_term = calc_logistic_vul(vul_par, n_age, prior);
  } else if(vul_type_term == "free") {
    for(int a=0;a<n_age-1;a++) vul_term(a) = exp(vul_par(a));
    vul_term(n_age-1) = Fratio * vul_term(n_age-2);
  } else {
    vul_term = calc_dome_vul(vul_par, n_age, prior);
  }

  // Terminal Year
  for(int a=0;a<n_age;a++) {
    F_at_age(n_y-1,a) = vul_term(a) * Fterm;
    Z_at_age(n_y-1,a) = F_at_age(n_y-1,a) + M(a);
    N(n_y-1,a) = Z_at_age(n_y-1,a) * CAA_hist(n_y-1,a);
    N(n_y-1,a) /= 1 - exp(-Z_at_age(n_y-1,a));
    N(n_y-1,a) /= F_at_age(n_y-1,a);
  }
  
  // Backwards recursion of N
  for(int y=n_y-1;y>0;y--) {
    for(int a=1;a<n_age;a++) {
      if(a==n_age-1) {
        F_at_age(y-1,a-1) = Newton_VPA_F_plus(F_at_age(y,a), Fratio, M(a-1), M(a), CAA_hist(y-1,a-1), CAA_hist(y-1,a), N(y,a), n_itF);
      } else {
        F_at_age(y-1,a-1) = CppAD::CondExpGt(CAA_hist(y-1,a-1), Type(1e-4), 
                 Newton_VPA_F(F_at_age(y,a), M(a-1), CAA_hist(y-1,a-1), N(y,a), n_itF), 
                 Type(1e-4));
      }
      Z_at_age(y-1,a-1) = F_at_age(y-1,a-1) + M(a-1);

      N(y-1,a-1) = CppAD::CondExpGt(CAA_hist(y-1,a-1), Type(1e-4), 
        Z_at_age(y-1,a-1) * CAA_hist(y-1,a-1)/(1 - exp(-Z_at_age(y-1,a-1)))/F_at_age(y-1,a-1), 
        N(y,a) * exp(Z_at_age(y-1,a-1)));

      if(a==n_age-1) {
        F_at_age(y-1,a) = Fratio * F_at_age(y-1,a-1);
        Z_at_age(y-1,a) = F_at_age(y-1,a) + M(a);
        
        N(y-1,a) = Z_at_age(y-1,a) * CAA_hist(y-1,a);
        N(y-1,a) /= 1 - exp(-Z_at_age(y-1,a));
        N(y-1,a) /= F_at_age(y-1,a);
      }
    }
  }
  
  for(int y=0;y<n_y;y++) {
    vector<Type> Fvec = F_at_age.row(y);
    F(y) = max(Fvec);
    for(int a=0;a<n_age;a++) {
      vul(y,a) = F_at_age(y,a)/F(y);
      CAApred(y,a) = N(y,a) * F_at_age(y,a) * (1 - exp(-Z_at_age(y,a)))/Z_at_age(y,a);
      B(y) += N(y,a) * weight(a);
      VB(y) += N(y,a) * weight(a) * vul(y,a);
    }
  }

  // Calculate nuisance parameters and likelihood
  vector<Type> q(nsurvey);
  array<Type> IAA(n_y,n_age,nsurvey);
  matrix<Type> IN(n_y,nsurvey);
  matrix<Type> Itot(n_y,nsurvey);
  matrix<Type> Ipred(n_y,nsurvey);
  IN.setZero();
  Itot.setZero();
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      for(int a=0;a<n_age;a++) {
        if(I_vul.col(sur).sum() > 0) {
          IAA(y,a,sur) = I_vul(a,sur) * N(y,a);
        } else {
          IAA(y,a,sur) = vul(a) * N(y,a);
        }
        IN(y,sur) += IAA(y,a,sur);
        if(I_units(sur)) Itot(y,sur) += IAA(y,a,sur) * weight(a); // Biomass vulnerable to survey
      }
    }
    if(!I_units(sur)) Itot.col(sur) = IN.col(sur); // Abundance vulnerable to survey
    q(sur) = calc_q(I_hist, Itot, sur, sur, Ipred, abs_I, n_y); // This function updates Ipred
  }

  vector<Type> nll_comp(nsurvey);
  nll_comp.setZero();
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      if(LWT(sur) > 0 && !R_IsNA(asDouble(I_hist(y,sur))) && I_hist(y,sur) > 0) {
        nll_comp(sur) -= dnorm_(log(I_hist(y,sur)), log(Ipred(y,sur)), I_sd(y,sur), true);
      }
      if(y>n_y-n_vulpen) {
        for(int a=0;a<n_age;a++) {
          prior -= dnorm_(log(vul(y,a)), log(vul(y-1,a)), vulpen, true);
        }
      }
      if(y>n_y-n_Rpen) prior -= dnorm_(log(N(y,1)), log(N(y-1,1)), Rpen, true);
    }
    nll_comp(sur) *= LWT(sur);
  }

  Type nll = nll_comp.sum() + penalty + prior;

  ADREPORT(q);
  if(CppAD::Variable(log_Fterm)) ADREPORT(Fterm);
  if(CppAD::Variable(log_Fratio)) ADREPORT(Fratio);

  REPORT(vul_par);
  REPORT(vul_term);
  REPORT(vul);
  REPORT(F);
  REPORT(F_at_age);
  REPORT(N);
  REPORT(CAApred);
  REPORT(Ipred);
  REPORT(VB);
  REPORT(B);

  REPORT(q);
  REPORT(Fterm);
  REPORT(Fratio);

  REPORT(nll_comp);
  REPORT(penalty);
  REPORT(prior);
  REPORT(nll);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif

