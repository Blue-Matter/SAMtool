
#ifndef DD_hpp
#define DD_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type DD(objective_function<Type> *obj) {

  DATA_SCALAR(S0);
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
  DATA_STRING(SR_type);
  DATA_STRING(condition);
  DATA_VECTOR(I_lambda);
  DATA_INTEGER(nsurvey);
  DATA_INTEGER(fix_sigma);
  DATA_INTEGER(state_space);

  PARAMETER(R0x);
  PARAMETER(transformed_h);
  PARAMETER(log_q_effort);
  PARAMETER(U_equilibrium);
  PARAMETER(log_omega);
  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_rec_dev);

  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else h = exp(transformed_h);
  h += 0.2;
  Type R0 = exp(R0x)/rescale;
  Type q_effort = exp(log_q_effort);
  Type omega = exp(log_omega);
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);


  //--DECLARING DERIVED VALUES
  Type Spr0 = (S0 * Alpha/(1 - S0) + wk)/(1 - Rho * S0);
  Type B0 = R0 * Spr0;
  Type N0 = R0/(1 - S0);

  Type Arec;
  Type Brec;

  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Arec /= Spr0;
    Brec = 5*h - 1;
    Brec /= (1-h) * B0;
  } else {
    Arec = pow(5*h, 1.25);
    Arec /= Spr0;
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= B0;
  }

  //--DECLARING STORAGE VECTORS
  int ny_p = ny + 1;
  int ny_k = ny + k;
  vector<Type> B(ny_p);
  vector<Type> N(ny_p);
  vector<Type> R(ny_k);
  vector<Type> Rec_dev(ny - k);

  vector<Type> Surv(ny);
  vector<Type> Cpred(ny);
  matrix<Type> Ipred(ny,nsurvey);
  vector<Type> Sp(ny);
  vector<Type> U(ny);

  //--INITIALIZE
  Type Seq = S0 * (1 - U_equilibrium);
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

  Type Ceqpred = B(0) * U_equilibrium;

  Type penalty = 0; // Penalty to likelihood for high U > 0.95
  Type prior = -dnorm_(log(B(0)/B0), log(dep), Type(0.01), true); // Penalty for initial depletion to get the corresponding U_equilibrium

  for(int tt=0; tt<ny; tt++){
    if(condition == "catch") {
      U(tt) = CppAD::CondExpLt(1 - C_hist(tt)/B(tt), Type(0.025),
        1 - posfun(1 - C_hist(tt)/B(tt), Type(0.025), penalty), C_hist(tt)/B(tt));
    } else {
      U(tt) = CppAD::CondExpLt(exp(-q_effort * E_hist(tt)), Type(0.025),
        1 - posfun(exp(-q_effort * E_hist(tt)), Type(0.025), penalty), 1 - exp(-q_effort * E_hist(tt)));
    }
    Surv(tt) = S0 * (1 - U(tt));
    Cpred(tt) = U(tt) * B(tt);
    Sp(tt) = B(tt) - Cpred(tt);

    if(SR_type == "BH") {
      R(tt+k) = BH_SR(B(tt), h, R0, B0);
    } else {
      R(tt+k) = Ricker_SR(B(tt), h, R0, B0);
    }
    if(state_space && tt + k < ny) {
      Rec_dev(tt) = exp(log_rec_dev(tt) - 0.5 * tau * tau);
      R(tt + k) *= Rec_dev(tt);
    }

    B(tt+1) = Surv(tt) * (Alpha * N(tt) + Rho * B(tt)) + wk * R(tt+1);
    N(tt+1) = Surv(tt) * N(tt) + R(tt+1);
  }

  //--ARGUMENTS FOR NLL
  // Objective function
  //creates storage for nll and sets value to 0
  vector<Type> q(nsurvey);
  if(condition == "catch") q = calc_q(I_hist, B, N, Ipred, nsurvey, I_units);

  vector<Type> nll_comp(nsurvey + 1);
  nll_comp.setZero();

  for(int tt=0; tt<ny; tt++){
    if(condition == "effort") {
      if(C_hist(tt) > 0) nll_comp(0) -= dnorm(log(C_hist(tt)), log(Cpred(tt)), omega, true);
    } else {
      for(int sur=0;sur<nsurvey;sur++) {
        for(int tt=0;tt<ny;tt++) {
          if(I_lambda(sur) > 0 && !R_IsNA(asDouble(I_hist(tt,sur)))) {
            if(fix_sigma) {
              nll_comp(sur) -= dnorm_(log(I_hist(tt,sur)), log(Ipred(tt,sur)), I_sd(tt,sur), true);
            } else {
              nll_comp(sur) -= dnorm(log(I_hist(tt,sur)), log(Ipred(tt,sur)), sigma, true);
            }
          }
        }
        nll_comp(sur) *= I_lambda(sur);
      }

    }
    if(state_space && tt + k < ny) nll_comp(nsurvey) -= dnorm(log_rec_dev(tt), Type(0), tau, true);
  }

  //Summing individual nll and penalties
  Type nll = nll_comp.sum() + penalty + prior;

  //-------REPORTING-------//
  ADREPORT(R0);
  ADREPORT(h);
  if(condition == "effort") ADREPORT(q_effort);
  if(condition == "catch") ADREPORT(q);
  if(CppAD::Variable(log_omega)) ADREPORT(omega);
  if(CppAD::Variable(log_sigma)) ADREPORT(sigma);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(condition == "effort") REPORT(omega);
  if(!fix_sigma) REPORT(sigma);
  if(state_space) REPORT(tau);
  REPORT(nll_comp);
  REPORT(nll);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(Spr0);
  REPORT(Cpred);
  REPORT(Ceqpred);
  REPORT(Ipred);
  REPORT(q);
  if(condition == "effort") REPORT(q_effort);
  REPORT(B);
  REPORT(N);
  REPORT(R);
  REPORT(log_rec_dev);
  REPORT(Rec_dev);
  REPORT(U);
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
