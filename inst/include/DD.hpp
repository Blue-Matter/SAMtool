
#ifndef DD_hpp
#define DD_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type DD(objective_function<Type> *obj) {
  using namespace ns_DD;

  DATA_SCALAR(Alpha);
  DATA_SCALAR(Rho);
  DATA_INTEGER(ny);
  DATA_INTEGER(k);
  DATA_SCALAR(wk);
  DATA_VECTOR(C_hist);
  DATA_SCALAR(dep);
  DATA_SCALAR(rescale);
  DATA_MATRIX(I_hist);
  DATA_IVECTOR(I_units);
  DATA_MATRIX(I_sd);
  DATA_VECTOR(E_hist);
  DATA_VECTOR(MW_hist);
  DATA_STRING(SR_type);
  DATA_STRING(condition);
  DATA_VECTOR(LWT);
  DATA_INTEGER(nsurvey);
  DATA_INTEGER(fix_sigma);
  DATA_INTEGER(nit_F);
  DATA_INTEGER(state_space);
  DATA_IVECTOR(use_prior); // Boolean vector, whether to set a prior for R0, h, M, q (length of 3 + nsurvey)
  DATA_MATRIX(prior_dist); // Distribution of priors for R0, h, M, q (rows), columns indicate parameters of distribution calculated in R (see make_prior fn)

  PARAMETER(R0x);
  PARAMETER(transformed_h);
  PARAMETER(log_M);
  PARAMETER(log_q_effort);
  PARAMETER(F_equilibrium);
  PARAMETER(log_omega);
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
  Type q_effort = exp(log_q_effort);
  Type omega = exp(log_omega);
  Type sigma = exp(log_sigma);
  Type sigma_W = exp(log_sigma_W);
  Type tau = exp(log_tau);

  //--DECLARING DERIVED VALUES
  Type S0 = exp(-M);
  Type Spr0 = (S0 * Alpha/(1 - S0) + wk)/(1 - Rho * S0);
  Type B0 = R0 * Spr0;
  Type N0 = R0/(1 - S0);

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
  Type Arec = CR/Spr0;

  //--DECLARING STORAGE VECTORS
  int ny_p = ny + 1;
  int ny_k = ny + k;
  vector<Type> B(ny_p);
  vector<Type> N(ny_p);
  vector<Type> R(ny_k);
  vector<Type> Rec_dev(ny);

  vector<Type> Surv(ny);
  vector<Type> Cpred(ny);
  matrix<Type> Ipred(ny,nsurvey);
  vector<Type> MWpred(ny);
  //vector<Type> Sp(ny);
  vector<Type> F(ny);

  //--INITIALIZE
  Type Z_equilibrium = F_equilibrium + M;
  Type Seq = S0 * exp(-F_equilibrium);
  Type SprEq = (Seq * Alpha/(1 - Seq) + wk)/(1 - Rho * Seq);
  Type Req;
  if(SR_type == "BH") {
    Req = Arec * SprEq - 1;
  } else {
    Req = log(Arec * SprEq);
  }
  Req /= Brec * SprEq;

  B(0) = Req * SprEq;
  N(0) = Req/(1 - Seq);
  for(int tt=0;tt<k;tt++) R(tt) = Req;
  if(state_space) {
    Rec_dev(0) = exp(log_rec_dev(0) - 0.5 * tau * tau);
    R(0) *= Rec_dev(0);
  }
  
  Type Ceqpred = F_equilibrium * B(0) * (1 - Seq)/(F_equilibrium + M);

  Type penalty = 0; // Penalty to likelihood for high F > 3
  Type prior = -dnorm_(log(B(0)/B0), log(dep), Type(0.01), true); // Penalty for initial depletion to get the corresponding U_equilibrium

  for(int tt=0; tt<ny; tt++){
    if(condition == "catch") {
      F(tt) = Newton_F(C_hist, M, B, Type(3), tt, nit_F, penalty);
    } else {
      Type tmp = Type(3) - q_effort * E_hist(tt);
      F(tt) = CppAD::CondExpGt(tmp, Type(0), Type(3) - posfun(tmp, Type(0), penalty), q_effort * E_hist(tt));
    }
    Surv(tt) = S0 * exp(-F(tt));
    Cpred(tt) = F(tt) * B(tt) * (1 - Surv(tt))/(F(tt) + M);
    MWpred(tt) = B(tt)/N(tt);

    if(SR_type == "BH") {
      R(tt+k) = BH_SR(B(tt), h, R0, B0);
    } else {
      R(tt+k) = Ricker_SR(B(tt), h, R0, B0);
    }
    if(state_space && tt<ny-1) {
      Rec_dev(tt+1) = exp(log_rec_dev(tt+1) - 0.5 * tau * tau);
      R(tt+1) *= Rec_dev(tt+1);
    }

    B(tt+1) = Surv(tt) * (Alpha * N(tt) + Rho * B(tt)) + wk * R(tt+1);
    N(tt+1) = Surv(tt) * N(tt) + R(tt+1);
  }

  //--ARGUMENTS FOR NLL
  vector<Type> q(nsurvey);
  if(condition == "catch") q = calc_q(I_hist, B, N, Ipred, nsurvey, I_units, ny);

  vector<Type> nll_comp(nsurvey + 2);
  nll_comp.setZero();

  for(int tt=0; tt<ny; tt++){
    if(condition == "effort") {
      if(C_hist(tt) > 0) nll_comp(0) -= dnorm(log(C_hist(tt)), log(Cpred(tt)), omega, true);
    } else {
      for(int sur=0;sur<nsurvey;sur++) {
        for(int tt=0;tt<ny;tt++) {
          if(LWT(sur) > 0 && !R_IsNA(asDouble(I_hist(tt,sur)))) {
            if(fix_sigma) {
              nll_comp(sur) -= dnorm_(log(I_hist(tt,sur)), log(Ipred(tt,sur)), I_sd(tt,sur), true);
            } else {
              nll_comp(sur) -= dnorm(log(I_hist(tt,sur)), log(Ipred(tt,sur)), sigma, true);
            }
          }
        }
        nll_comp(sur) *= LWT(sur);
      }
    }
    if(LWT(nsurvey) > 0 && !R_IsNA(asDouble(MW_hist(tt)))) {
      nll_comp(nsurvey) -= LWT(nsurvey) * dnorm(log(MW_hist(tt)), log(MWpred(tt)), sigma_W, true);
    }
    if(state_space) nll_comp(nsurvey+1) -= dnorm(log_rec_dev(tt), Type(0), tau, true);
  }

  //Summing individual nll and penalties
  prior -= calc_prior(use_prior, prior_dist, R0, h, SR_type == "BH", log_M, q);
  Type nll = nll_comp.sum() + penalty + prior;

  //-------REPORTING-------//
  ADREPORT(R0);
  ADREPORT(h);
  if(CppAD::Variable(log_M)) ADREPORT(M);
  if(condition == "effort") ADREPORT(q_effort);
  if(condition == "catch") ADREPORT(q);
  if(CppAD::Variable(log_omega)) ADREPORT(omega);
  if(CppAD::Variable(log_sigma)) ADREPORT(sigma);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(MW_hist.sum() > 0) {
    ADREPORT(sigma_W);
    REPORT(sigma_W);
  }
  if(condition == "effort") REPORT(omega);
  if(!fix_sigma) REPORT(sigma);
  if(state_space) REPORT(tau);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(CR);
  REPORT(Spr0);
  REPORT(Cpred);
  REPORT(Ceqpred);
  REPORT(Ipred);
  REPORT(q);
  if(condition == "effort") REPORT(q_effort);
  REPORT(B);
  REPORT(N);
  REPORT(MWpred);
  REPORT(R);
  REPORT(log_rec_dev);
  REPORT(Rec_dev);
  REPORT(F);
  REPORT(M);
  REPORT(h);
  REPORT(R0);
  REPORT(N0);
  REPORT(B0);
  REPORT(penalty);
  REPORT(prior);

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
