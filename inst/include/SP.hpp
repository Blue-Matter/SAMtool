
#ifndef SP_hpp
#define SP_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type SP(objective_function<Type> *obj) {

  using namespace ns_SP;

  DATA_VECTOR(C_hist);
  DATA_SCALAR(rescale);
  DATA_MATRIX(I_hist);
  DATA_MATRIX(I_sd);
  DATA_VECTOR(I_lambda);
  DATA_INTEGER(fix_sigma);
  DATA_INTEGER(nsurvey);
  DATA_INTEGER(ny);
  DATA_IVECTOR(est_B_dev);
  DATA_INTEGER(nstep);
  DATA_SCALAR(dt);
  DATA_INTEGER(n_itF);
  DATA_IVECTOR(use_prior); // Boolean vector, whether to set a prior for r, MSY
  DATA_MATRIX(prior_dist); // Distribution of priors, columns indicate parameters of distribution calculated in R (see make_prior_SP fn)
  DATA_INTEGER(sim_process_error);
  //DATA_VECTOR_INDICATOR(keep, I_hist);

  PARAMETER(log_FMSY);
  PARAMETER(MSYx);
  PARAMETER(log_dep);
  PARAMETER(log_n);
  PARAMETER_VECTOR(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_B_dev);

  Type FMSY = exp(log_FMSY);
  Type MSY = exp(MSYx)/rescale;
  Type dep = exp(log_dep);
  Type n = exp(log_n);

  vector<Type> sigma(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) sigma(sur) = exp(log_sigma(sur));
  Type tau = exp(log_tau);

  Type euler = exp(Type(1.0));

  Type n_term = CppAD::CondExpEq(n, Type(1), euler, pow(n, n/(n-1)));
  Type n_term2 = CppAD::CondExpEq(n, Type(1), Type(1/euler), pow(n, 1/(1-n)));

  Type BMSY = MSY/FMSY;
  Type K = BMSY / n_term2;
  Type r = MSY * n_term / K; // r = FMSY * n

  vector<Type> B(ny+1);
  vector<Type> B_dev(ny);
  vector<Type> SP(ny);
  vector<Type> Cpred(ny);
  matrix<Type> Ipred(ny,nsurvey);
  vector<Type> F(ny);
  vector<Type> log_B_dev_sim = log_B_dev;
  
  B_dev.fill(1);

  Type penalty = 0;
  Type prior = 0;

  if(use_prior(0)) { // log-normal r prior with log-Jacobian transformation = 0, exact with fixed n
    prior -= dnorm_(log(r), prior_dist(0,0), prior_dist(0,1), true);
  }
  if(use_prior(1)) { // log-normal MSY prior with log-Jacobian transformation = 0
    prior -= dnorm_(log(MSY), prior_dist(1,0), prior_dist(1,1), true);
  }

  B(0) = dep * K;
  for(int y=0;y<ny;y++) {
    if(est_B_dev(y)) {
      B_dev(y) = exp(log_B_dev(y) - 0.5 * tau * tau);
      SIMULATE if(sim_process_error) {
        log_B_dev_sim(y) = rnorm(log_B_dev(y), tau);
        B_dev(y) = exp(log_B_dev_sim(y) - 0.5 * tau * tau);
      }
      B(y) *= B_dev(y);
    }
    
    if(C_hist(y) > 1e-8) {
      F(y) = SP_F(C_hist(y)/(C_hist(y) + B(y)), C_hist(y), MSY, K, n, n_term, dt, nstep, n_itF, Cpred, B, y, penalty);
    } else {
      F(y) = SP_F(C_hist(y)/(C_hist(y) + B(y)), C_hist(y), MSY, K, n, n_term, dt, 1, n_itF, Cpred, B, y, penalty);
    }
    SP(y) = B(y+1) - B(y) + Cpred(y);
    
    SIMULATE {
      C_hist(y) = Cpred(y);
    }
  }

  vector<Type> q = calc_q(I_hist, B, Ipred, nsurvey);

  vector<Type> nll_comp(I_hist.cols()+1);
  nll_comp.setZero();

  for(int sur=0;sur<I_hist.cols();sur++) {
    for(int y=0;y<ny;y++) {
      if(I_lambda(sur) > 0 && !R_IsNA(asDouble(I_hist(y,sur)))) {
        if(fix_sigma) {
          nll_comp(sur) -= dnorm_(log(I_hist(y,sur)), log(Ipred(y,sur)), I_sd(y,sur), true);
          SIMULATE {
            I_hist(y,sur) = exp(rnorm(log(Ipred(y,sur)), I_sd(y,sur)));
          }
        } else {
          nll_comp(sur) -= dnorm(log(I_hist(y,sur)), log(Ipred(y,sur)), sigma(sur), true);
          SIMULATE {
            I_hist(y,sur) = exp(rnorm(log(Ipred(y,sur)), sigma(sur)));
          }
        }
      }
    }
    nll_comp(sur) *= I_lambda(sur);
  }
  for(int y=0;y<ny;y++) {
    if(est_B_dev(y)) nll_comp(nsurvey) -= dnorm(log_B_dev(y), Type(0), tau, true);
  }

  Type nll = nll_comp.sum() + penalty + prior;

  Type F_FMSY_final = F(F.size()-1)/FMSY;
  Type B_BMSY_final = B(B.size()-1)/BMSY;
  Type B_K_final = B(B.size()-1)/K;

  ADREPORT(FMSY);
  ADREPORT(MSY);
  if(CppAD::Variable(log_dep)) ADREPORT(dep);
  if(CppAD::Variable(log_n)) ADREPORT(n);
  ADREPORT(q);
  ADREPORT(r);
  ADREPORT(K);
  if(!fix_sigma) ADREPORT(sigma);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  ADREPORT(F_FMSY_final);
  ADREPORT(B_BMSY_final);
  ADREPORT(B_K_final);

  REPORT(FMSY);
  REPORT(MSY);
  REPORT(dep);
  REPORT(n);
  REPORT(n_term);
  REPORT(q);
  REPORT(sigma);
  REPORT(tau);
  REPORT(r);
  REPORT(K);
  REPORT(BMSY);
  REPORT(Cpred);
  REPORT(Ipred);
  REPORT(B);
  REPORT(SP);
  REPORT(F);
  REPORT(log_B_dev);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);
  
  SIMULATE {
    REPORT(C_hist);
    REPORT(I_hist);
    REPORT(log_B_dev_sim);
  }

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
