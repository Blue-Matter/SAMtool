
#ifndef cDD_hpp
#define cDD_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type cDD(objective_function<Type> *obj) {
  using namespace ns_cDD;

  DATA_SCALAR(M);
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
  DATA_STRING(SR_type);
  DATA_INTEGER(nitF);
  DATA_VECTOR(I_lambda);
  DATA_INTEGER(nsurvey);
  DATA_INTEGER(fix_sigma);
  DATA_INTEGER(state_space);

  PARAMETER(R0x);
  PARAMETER(transformed_h);
  PARAMETER(F_equilibrium);
  PARAMETER(log_sigma);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_rec_dev);

  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else h = exp(transformed_h);
  h += 0.2;
  Type R0 = exp(R0x)/rescale;
  Type sigma = exp(log_sigma);
  Type tau = exp(log_tau);
  int SR_type2 = SR_type == "BH";


  //--DECLARING DERIVED VALUES
  Type BPR0 = cDD_BPR(Type(0), M, wk, Kappa, Winf);
  Type B0 = BPR0 * R0;
  Type N0 = R0/M;

  Type Arec;
  Type Brec;

  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Arec /= BPR0;
    Brec = 5*h - 1;
    Brec /= (1-h) * B0;
  } else {
    Arec = pow(5*h, 1.25);
    Arec /= BPR0;
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= B0;
  }

  //--DECLARING STORAGE VECTORS
  vector<Type> B(ny+1);
  vector<Type> N(ny+1);
  vector<Type> R(ny+k);

  vector<Type> F(ny);
  vector<Type> Z(ny);
  vector<Type> Cpred(ny);
  matrix<Type> Ipred(ny,nsurvey);

  vector<Type> BPRinf(ny);
  vector<Type> Rec_dev(ny-k);
  vector<Type> Binf(ny);
  vector<Type> Ninf(ny);

  //--INITIALIZE
  Type BPReq = cDD_BPR(F_equilibrium, M, wk, Kappa, Winf);
  Type Req = cDD_R(BPReq, Arec, Brec, SR_type2);

  B(0) = Req * BPReq;
  N(0) = Req/(F_equilibrium + M);
  for(int tt=0;tt<k;tt++) R(tt) = Req;

  Type Ceqpred = F_equilibrium * B(0);

  Type penalty = 0; // Pentalty to likelihood for high F > max_F
  Type prior = -dnorm_(log(B(0)/B0), log(dep), Type(0.01), true); // Penalty for initial depletion to get corresponding F

  for(int tt=0; tt<ny; tt++) {
    Type F_start = CppAD::CondExpLe(C_hist(tt), Type(1e-8), Type(0), -log(1 - C_hist(tt)/B(tt)));
    F(tt) = cDD_F(F_start, C_hist(tt), M, Winf, Kappa, wk, N, B, Cpred, BPRinf, Binf, R, Ninf,
      CppAD::Integer(CppAD::CondExpLe(C_hist(tt), Type(1e-8), Type(1), Type(nitF))), tt);
    Z(tt) = F(tt) + M;

    N(tt+1) = Ninf(tt) + (N(tt) - Ninf(tt)) * exp(-Z(tt));
    B(tt+1) = Binf(tt);
    B(tt+1) += Kappa * Winf * (N(tt) - Ninf(tt)) / (Z(tt) + Kappa);
    B(tt+1) += (B(tt) - Binf(tt) - Kappa * Winf * (N(tt) - Ninf(tt))/(Z(tt) + Kappa)) * exp(-Z(tt) - Kappa);

    if(SR_type == "BH") {
      R(tt+k) = BH_SR(B(tt), h, R0, B0);
    } else {
      R(tt+k) = Ricker_SR(B(tt), h, R0, B0);
    }
    if(state_space && tt+k<ny) {
      Rec_dev(tt) = exp(log_rec_dev(tt) - 0.5 * tau * tau);
      R(tt+k) *= Rec_dev(tt);
    }
  }

  //--ARGUMENTS FOR NLL
  vector<Type> q = calc_q(I_hist, B, N, Ipred, nsurvey, I_units);

  // Objective function
  //creates storage for jnll and sets value to 0
  vector<Type> nll_comp(nsurvey+1);
  nll_comp.setZero();

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
  if(state_space) {
    for(int tt=0;tt<log_rec_dev.size();tt++) nll_comp(nsurvey) -= dnorm(log_rec_dev(tt), Type(0), tau, true);
  }

  //Summing individual jnll and penalties
  Type nll = nll_comp.sum() + penalty + prior;

  //-------REPORTING-------//
  ADREPORT(R0);
  ADREPORT(h);
  ADREPORT(q);
  if(CppAD::Variable(log_sigma)) ADREPORT(sigma);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(!fix_sigma) REPORT(sigma);
  if(state_space) REPORT(tau);
  REPORT(nll);
  REPORT(Arec);
  REPORT(Brec);
  REPORT(BPR0);
  REPORT(Cpred);
  REPORT(Ceqpred);
  REPORT(Ipred);
  REPORT(q);
  REPORT(B);
  REPORT(N);
  REPORT(R);
  REPORT(F);
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

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
