
#ifndef SRA_scope_hpp
#define SRA_scope_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type SRA_scope(objective_function<Type> *obj) {

  using namespace ns_SRA_scope;

  DATA_MATRIX(C_hist);    // Total catch by year and fleet
  DATA_VECTOR(C_eq);      // Equilibrium catch by fleet

  DATA_MATRIX(E_hist);    // Effort by year and fleet
  DATA_VECTOR(E_eq);      // Equilibrium effort by fleet

  DATA_STRING(condition); // Indicates whether the model will condition on 'effort', 'catch' (estimated F's), 'catch2' (internally optimize for F)
  DATA_INTEGER(nll_C);    // Indicates whether there is a likelihood for the catch. TRUE when (a) condition = 'catch' or (b) 'effort' with nfleet > 1

  DATA_MATRIX(I_hist);    // Index by year and survey
  DATA_MATRIX(sigma_I);   // Standard deviation of index by year and survey

  DATA_ARRAY(CAA_hist);   // Catch-at-age re-weighted by year, age, fleet
  DATA_MATRIX(CAA_n);     // Annual samples in CAA by year and fleet

  DATA_ARRAY(CAL_hist);   // Catch-at-length re-weighted by year, length_bin, fleet
  DATA_MATRIX(CAL_n);     // Annual samples in CAL by year and fleet

  DATA_ARRAY(s_CAA_hist);  // Catch-at-age re-weighted by year, age, survey
  DATA_MATRIX(s_CAA_n);    // Annual samples in CAA by year and survey

  DATA_ARRAY(s_CAL_hist);  // Catch-at-length re-weighted by year, length_bin, survey
  DATA_MATRIX(s_CAL_n);    // Annual samples in CAL by year and survey

  DATA_VECTOR(length_bin); // Vector of length bins
  DATA_MATRIX(msize);      // Vector of annual mean size by year and fleet
  DATA_STRING(msize_type); // Whether the mean size is length or weight

  DATA_IMATRIX(sel_block); // Assigns selectivity by year and fleet
  DATA_INTEGER(nsel_block); // The number of selectivity "blocks"

  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(max_age);  // Maximum age (plus-group)
  DATA_INTEGER(nfleet);   // Number of fleets
  DATA_INTEGER(nsurvey);  // Number of surveys

  DATA_MATRIX(M);         // Natural mortality at age
  DATA_MATRIX(len_age);   // Length-at-age
  DATA_SCALAR(Linf);      // Linf
  DATA_SCALAR(CV_LAA);    // CV of length-at-age
  DATA_MATRIX(wt);        // Weight-at-age
  DATA_MATRIX(mat);       // Maturity-at-age at the beginning of the year

  DATA_IVECTOR(vul_type); // Integer vector indicating whether free (-2), logistic (-1), or dome vul (0) is used
  DATA_IVECTOR(s_vul_type); // Same but for surveys, but can mirror to B (-4), SSB (-3), or fleet (>0)
  DATA_IVECTOR(abs_I);    // Boolean, whether index is an absolute (fix q = 1) or relative terms (estimate q)
  DATA_IVECTOR(I_units);  // Boolean, whether index is biomass based (= 1) or abundance-based (0)

  DATA_MATRIX(age_error); // Ageing error matrix

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_MATRIX(LWT_C);     // LIkelihood weights for catch, CAA, CAL, MS, C_eq
  DATA_MATRIX(LWT_Index); // Likelihood weights for the index
  DATA_STRING(comp_like); // Whether to use "multinomial" or "lognormal" distribution for age/lengthc comps

  DATA_SCALAR(max_F);     // Maximum F in the model
  DATA_SCALAR(rescale);   // R0 rescaler
  DATA_INTEGER(ageM);     // Age of maturity used for averaging E0 and EPR0

  DATA_IVECTOR(est_early_rec_dev); // Indicates years in which log_early_rec_dev are estimated. Then, lognormal bias correction estimates are added..
  DATA_IVECTOR(est_rec_dev); // Indicates years in which log_rec_dev are estimated.
  DATA_IVECTOR(yind_F);   // When condition = "catch", the year in F's are estimated and all other F parameters are deviations from this F
  DATA_INTEGER(nit_F);    // When condition = "catch2", the number of iterations for Newton-Raphson method to solve for F
  DATA_INTEGER(plusgroup) // Boolean, whether the maximum age in the plusgroup is modeled.

  PARAMETER(R0x);                       // Unfished recruitment
  PARAMETER(transformed_h);             // Steepness
  PARAMETER_MATRIX(vul_par);            // Matrix of vul_par 3 rows and nsel_block columns
  PARAMETER_MATRIX(s_vul_par);          // Matrix of selectivity parameters, 3 rows and nsurvey columns
  PARAMETER_VECTOR(log_q_effort);       // log_q for F when condition = "effort"
  PARAMETER_MATRIX(log_F);              // log_F_deviations when condition = "catch"
  PARAMETER_VECTOR(log_F_equilibrium);  // Equilibrium F by fleet when condition != "effort"

  PARAMETER_VECTOR(log_CV_msize);       // CV of mean size
  PARAMETER(log_tau);                   // Std. dev. of rec devs
  PARAMETER_VECTOR(log_early_rec_dev);  // Rec devs for first year abundance
  PARAMETER_VECTOR(log_rec_dev);        // Rec devs for all other years in the model

  int nlbin = length_bin.size();
  Type bin_width = length_bin(1) - length_bin(0);

  Type R0 = exp(R0x)/rescale;
  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else {
    h = exp(transformed_h);
  }
  h += 0.2;

  Type tau = exp(log_tau);

  // Vulnerability (length-based) and F parameters
  Type penalty = 0;
  Type prior = 0.;
  vector<Type> LFS(nsel_block);
  vector<Type> L5(nsel_block);
  vector<Type> Vmaxlen(nsel_block);
  array<Type> vul = calc_vul(vul_par, vul_type, len_age, LFS, L5, Vmaxlen, Linf, nfleet, sel_block, nsel_block, prior);

  vector<Type> q_effort(nfleet);
  vector<Type> CV_msize(nfleet);
  vector<Type> F_equilibrium(nfleet);
  F_equilibrium.setZero();

  matrix<Type> F(n_y,nfleet);
  matrix<Type> Z = M;

  for(int ff=0;ff<nfleet;ff++) {
    CV_msize(ff) = exp(log_CV_msize(ff));
    q_effort(ff) = exp(log_q_effort(ff));
    if(condition != "effort" && C_eq(ff)>0) F_equilibrium(ff) = exp(log_F_equilibrium(ff));
    if(condition == "effort" && E_eq(ff)>0) F_equilibrium(ff) = q_effort(ff) * E_eq(ff);
    if(condition == "catch") {
      Type tmp = max_F - exp(log_F(yind_F(ff),ff));
      F(yind_F(ff),ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), exp(log_F(yind_F(ff),ff)));
    }
  }

  ////// Equilibrium reference points and per-recruit quantities - calculate annually
  vector<vector<Type> > NPR_unfished(n_y);
  vector<Type> EPR0(n_y);
  vector<Type> E0(n_y);
  vector<Type> B0(n_y);
  vector<Type> N0(n_y);

  Type E0_SR = 0;
  for(int y=0;y<n_y;y++) {
    NPR_unfished(y) = calc_NPR0(M, max_age, y, plusgroup);

    EPR0(y) = sum_EPR(NPR_unfished(y), wt, mat, max_age, y);
    E0(y) = R0 * EPR0(y);
    B0(y) = R0 * sum_BPR(NPR_unfished(y), wt, max_age, y);
    N0(y) = R0 * NPR_unfished(y).sum();

    if(y < ageM) E0_SR += E0(y);
  }
  E0_SR /= Type(ageM);
  Type EPR0_SR = E0_SR/R0;

  Type Arec, Brec;
  if(SR_type == "BH") {
    Arec = 4 *h;
    Arec /= 1-h;
    Brec = 5*h - 1;
    Brec /= (1-h);
  } else {
    Arec = pow(5*h, 1.25);
    Brec = 1.25;
    Brec *= log(5*h);
  }
  Arec /= EPR0_SR;
  Brec /= E0_SR;

  ////// During time series year = 1, 2, ..., n_y
  vector<matrix<Type> > ALK(n_y);
  matrix<Type> N(n_y+1, max_age);

  vector<Type> C_eq_pred(nfleet);
  array<Type> CAAtrue(n_y, max_age, nfleet);   // Catch (in numbers) at year and age at the mid-point of the season
  array<Type> CAApred(n_y, max_age, nfleet);   // Catch (in numbers) at year and age at the mid-point of the season
  array<Type> CALpred(n_y, nlbin, nfleet);
  matrix<Type> MLpred(n_y, nfleet);
  matrix<Type> MWpred(n_y, nfleet);
  matrix<Type> CN(n_y, nfleet);             // Catch in numbers

  matrix<Type> Cpred(n_y, nfleet);
  matrix<Type> Ipred(n_y, nsurvey);          // Predicted index at year

  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(max_age-1);
  matrix<Type> VB(n_y+1, nfleet);   // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year

  C_eq_pred.setZero();
  CAApred.setZero();
  CALpred.setZero();
  MLpred.setZero();
  CN.setZero();

  Cpred.setZero();
  Ipred.setZero();

  VB.setZero();
  B.setZero();
  E.setZero();

  // Equilibrium quantities (leading into first year of model)
  vector<Type> NPR_equilibrium = calc_NPR(F_equilibrium, vul, nfleet, M, max_age, 0, plusgroup);
  Type EPR_eq = sum_EPR(NPR_equilibrium, wt, mat, max_age, 0);
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

    B(0) += N(0,a) * wt(0,a);
    E(0) += N(0,a) * wt(0,a) * mat(0,a);

    Type Z_eq = M(0,a);
    for(int ff=0;ff<nfleet;ff++) Z_eq += vul(0,a,ff) * F_equilibrium(ff);
    Type mean_N_eq = N(0,a) * (1 - exp(-Z_eq)) / Z_eq;
    for(int ff=0;ff<nfleet;ff++) {
      C_eq_pred(ff) += vul(0,a,ff) * F_equilibrium(ff) * mean_N_eq * wt(0,a);
      VB(0,ff) += N(0,a) * wt(0,a) * vul(0,a,ff);
    }
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y), h, R0, E0_SR);
    } else {
      R(y+1) = Ricker_SR(E(y), h, R0, E0_SR);
    }

    if(y<n_y-1 && est_rec_dev(y+1)) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * tau * tau);

    N(y+1,0) = R(y+1);
    if(Type(max_age) != Linf) ALK(y) = generate_ALK(length_bin, len_age, CV_LAA, max_age, nlbin, bin_width, y);

    if(condition == "catch") {
      for(int ff=0;ff<nfleet;ff++) {
        if(y != yind_F(ff)) {
          Type tmp = max_F - F(yind_F(ff),ff) * exp(log_F(y,ff));
          F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), F(yind_F(ff),ff) * exp(log_F(y,ff)));
        }
      }
    } else if(condition == "catch2") {
      F.row(y) = Newton_SRA_F(C_hist, N, M, wt, VB, vul, max_F, y, max_age, nfleet, nit_F, penalty);
    } else {
      for(int ff=0;ff<nfleet;ff++) {
        Type tmp = max_F - q_effort(ff) * E_hist(y,ff);
        F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), q_effort(ff) * E_hist(y,ff));
      }
    }

    for(int a=0;a<max_age;a++) {
      for(int ff=0;ff<nfleet;ff++) Z(y,a) += vul(y,a,ff) * F(y,ff);
      Type mean_N = N(y,a) * (1 - exp(-Z(y,a))) / Z(y,a);

      if(a<max_age-1) N(y+1,a+1) = N(y,a) * exp(-Z(y,a));
      if(plusgroup && a==max_age-1) N(y+1,a) += N(y,a) * exp(-Z(y,a));

      for(int ff=0;ff<nfleet;ff++) {
        CAAtrue(y,a,ff) = vul(y,a,ff) * F(y,ff) * mean_N;
        CN(y,ff) += CAAtrue(y,a,ff);
        Cpred(y,ff) += CAAtrue(y,a,ff) * wt(y,a);

        if(Type(max_age) != Linf) {
          for(int len=0;len<nlbin;len++) {
            CALpred(y,len,ff) += CAAtrue(y,a,ff) * ALK(y)(a,len);
            MLpred(y,ff) += CAAtrue(y,a,ff) * ALK(y)(a,len) * length_bin(len);
          }
        }
        for(int aa=0;aa<max_age;aa++) CAApred(y,aa,ff) += CAAtrue(y,a,ff) * age_error(a,aa); // a = true, aa = observed ages
        VB(y+1,ff) += vul(y+1,a,ff) * N(y+1,a) * wt(y+1,a);
      }

      B(y+1) += N(y+1,a) * wt(y+1,a);
      E(y+1) += N(y+1,a) * wt(y+1,a) * mat(y+1,a);
    }
    if(Type(max_age) != Linf) for(int ff=0;ff<nfleet;ff++) MLpred(y,ff) /= CN(y,ff);
    if(msize_type == "weight") for(int ff=0;ff<nfleet;ff++) MWpred(y,ff) = Cpred(y,ff)/CN(y,ff);
  }

  // Calculate nuisance parameters and likelihood
  // Survey selectivity
  vector<Type> s_LFS(nsurvey);
  vector<Type> s_L5(nsurvey);
  vector<Type> s_Vmaxlen(nsurvey);

  array<Type> s_CAAtrue(n_y, max_age, nsurvey); // True abundance at age vulnerable to survey
  array<Type> s_CAApred(n_y, max_age, nsurvey); // Predicted abundance (after ageing error) at age vulnerable to survey
  array<Type> s_CALpred(n_y, nlbin, nsurvey); // Abundance at length vulnerable to survey
  matrix<Type> s_CN(n_y, nsurvey); // Total abundance vulnerable to the survey
  matrix<Type> s_BN(n_y, nsurvey); // Biomass or abundance vulnerable to the survey

  s_CAApred.setZero();
  s_CALpred.setZero();
  s_CN.setZero();
  s_BN.setZero();

  array<Type> s_vul = calc_vul_sur(s_vul_par, s_vul_type, len_age, s_LFS, s_L5, s_Vmaxlen, Linf, mat, vul, prior);
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      for(int a=0;a<max_age;a++) {
        s_CAAtrue(y,a,sur) = s_vul(y,a,sur) * N(y,a);
        s_CN(y,sur) += s_CAAtrue(y,a,sur);

        for(int aa=0;aa<max_age;aa++) s_CAApred(y,aa,sur) += s_CAAtrue(y,a,sur) * age_error(a,aa);

        if(I_units(sur)) s_BN(y,sur) += s_CAAtrue(y,a,sur) * wt(y,a); // Biomass vulnerable to survey
        if(Type(max_age) != Linf && s_CAL_n.col(sur).sum() > 0) {
          for(int len=0;len<nlbin;len++) s_CALpred(y,len,sur) += s_CAAtrue(y,a,sur) * ALK(y)(a,len);
        }
      }
    }
    if(!I_units(sur)) s_BN.col(sur) = s_CN.col(sur); // Abundance vulnerable to survey
    q(sur) = calc_q(I_hist, s_BN, sur, sur, Ipred, abs_I);
  }

  vector<Type> nll_Catch(nfleet);
  vector<Type> nll_Index(nsurvey);
  vector<Type> nll_s_CAA(nsurvey);
  vector<Type> nll_s_CAL(nsurvey);
  vector<Type> nll_CAA(nfleet);
  vector<Type> nll_CAL(nfleet);
  vector<Type> nll_MS(nfleet);
  Type nll_log_rec_dev = 0;
  vector<Type> nll_Ceq(nfleet);

  nll_Catch.setZero();
  nll_Index.setZero();
  nll_CAA.setZero();
  nll_CAL.setZero();
  nll_s_CAA.setZero();
  nll_s_CAL.setZero();
  nll_MS.setZero();
  nll_Ceq.setZero();

  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      if(LWT_Index(sur,0) > 0 && !R_IsNA(asDouble(I_hist(y,sur)))) {
        nll_Index(sur) -= dnorm_(log(I_hist(y,sur)), log(Ipred(y,sur)), sigma_I(y,sur), true);
      }

      if(LWT_Index(sur,1) > 0 && !R_IsNA(asDouble(s_CAA_n(y,sur))) && s_CAA_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_s_CAA(sur) -= comp_multinom(s_CAA_hist, s_CAApred, s_CN, s_CAA_n, y, max_age, sur);
        } else {
          nll_s_CAA(sur) -= comp_lognorm(s_CAA_hist, s_CAApred, s_CN, s_CAA_n, y, max_age, sur);
        }
      }

      if(LWT_Index(sur,2) > 0 && !R_IsNA(asDouble(s_CAL_n(y,sur))) && s_CAL_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_s_CAL(sur) -= comp_multinom(s_CAL_hist, s_CALpred, s_CN, s_CAL_n, y, nlbin, sur);
        } else {
          nll_s_CAL(sur) -= comp_lognorm(s_CAL_hist, s_CALpred, s_CN, s_CAL_n, y, nlbin, sur);
        }
      }
    }
    nll_Index(sur) *= LWT_Index(sur,0);
    nll_s_CAA(sur) *= LWT_Index(sur,1);
    nll_s_CAL(sur) *= LWT_Index(sur,2);
  }

  for(int ff=0;ff<nfleet;ff++) {
    for(int y=0;y<n_y;y++) {
      if(C_hist(y,ff)>0 || E_hist(y,ff)>0) {
        if(LWT_C(ff,1) > 0 && !R_IsNA(asDouble(CAA_n(y,ff))) && CAA_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_CAA(ff) -= comp_multinom(CAA_hist, CAApred, CN, CAA_n, y, max_age, ff);
          } else {
            nll_CAA(ff) -= comp_lognorm(CAA_hist, CAApred, CN, CAA_n, y, max_age, ff);
          }
        }

        if(LWT_C(ff,2) > 0 && !R_IsNA(asDouble(CAL_n(y,ff))) && CAL_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_CAL(ff) -= comp_multinom(CAL_hist, CALpred, CN, CAL_n, y, nlbin, ff);
          } else {
            nll_CAL(ff) -= comp_lognorm(CAL_hist, CALpred, CN, CAL_n, y, nlbin, ff);
          }
        }

        if(nll_C && LWT_C(ff,0) > 0) nll_Catch(ff) -= dnorm_(log(C_hist(y,ff)), log(Cpred(y,ff)), Type(0.01), true);
        if(LWT_C(ff,3) > 0 && !R_IsNA(asDouble(msize(y,ff))) && msize(y,ff) > 0) {
          if(msize_type == "length") {
            nll_MS(ff) -= dnorm_(msize(y,ff), MLpred(y,ff), CV_msize(ff) * msize(y,ff), true);
          } else {
            nll_MS(ff) -= dnorm_(msize(y,ff), MWpred(y,ff), CV_msize(ff) * msize(y,ff), true);
          }
        }
      }
    }

    nll_Catch(ff) *= LWT_C(ff,0);
    nll_CAA(ff) *= LWT_C(ff,1);
    nll_CAL(ff) *= LWT_C(ff,2);
    nll_MS(ff) *= LWT_C(ff,3);

    if(LWT_C(ff,4) > 0 && C_eq(ff) > 0 && condition != "effort") {
      nll_Ceq(ff) = -1 * LWT_C(ff,4) * dnorm_(log(C_eq(ff)), log(C_eq_pred(ff)), Type(0.01), true);
    }
  }

  for(int y=0;y<n_y;y++) {
    if(est_rec_dev(y)) nll_log_rec_dev -= dnorm_(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<max_age-1;a++) {
    if(est_early_rec_dev(a)) nll_log_rec_dev -= dnorm_(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_Catch.sum() + nll_Index.sum();
  nll += nll_s_CAA.sum() + nll_s_CAL.sum();
  nll += nll_CAA.sum() + nll_CAL.sum() + nll_MS.sum();
  nll += nll_log_rec_dev + nll_Ceq.sum();
  nll += penalty + prior;

  if(CppAD::Variable(R0x)) ADREPORT(R0);
  ADREPORT(h);
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(condition == "effort") ADREPORT(q_effort);
  if(nll_Index.sum() != 0) ADREPORT(q);

  REPORT(R0x);
  REPORT(transformed_h);
  REPORT(vul_par);
  REPORT(LFS);
  REPORT(L5);
  REPORT(Vmaxlen);
  REPORT(log_q_effort);
  REPORT(log_F_equilibrium);

  if(nll_MS.sum() != 0) REPORT(log_CV_msize);
  REPORT(log_tau);
  REPORT(log_early_rec_dev);
  REPORT(log_rec_dev);

  REPORT(R0);
  REPORT(h);
  REPORT(tau);
  if(nll_MS.sum() != 0) REPORT(CV_msize);
  if(condition == "catch") REPORT(log_F);
  REPORT(F_equilibrium);
  REPORT(vul);
  REPORT(F);
  REPORT(Z);

  REPORT(NPR_unfished);
  REPORT(EPR0);
  REPORT(B0);
  REPORT(E0);
  REPORT(N0);

  REPORT(Arec);
  REPORT(Brec);
  REPORT(E0_SR);
  REPORT(EPR0_SR);

  if(nll_CAL.sum() != 0 || nll_s_CAL.sum() != 0 || ((nll_MS.sum() != 0) & (msize_type == "length"))) REPORT(ALK);
  REPORT(N);
  REPORT(CAApred);
  REPORT(CALpred);
  REPORT(MLpred);
  if(msize_type == "weight") REPORT(MWpred);
  REPORT(CN);
  REPORT(Cpred);
  REPORT(Ipred);
  REPORT(R);
  REPORT(R_early);
  REPORT(VB);
  REPORT(B);
  REPORT(E);

  REPORT(NPR_equilibrium);
  REPORT(EPR_eq);
  REPORT(R_eq);

  REPORT(C_eq_pred);
  if(nll_Index.sum() != 0) REPORT(q);

  REPORT(nll_Catch);
  REPORT(nll_Index);
  REPORT(nll_s_CAA);
  REPORT(nll_s_CAL);
  REPORT(nll_CAA);
  REPORT(nll_CAL);
  REPORT(nll_MS);
  REPORT(nll_Ceq);
  REPORT(nll_log_rec_dev);

  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  if(age_error.trace() != Type(max_age)) {
    REPORT(CAAtrue);
    REPORT(s_CAAtrue);
  }

  if(nll_Index.sum() != 0) {
    REPORT(s_vul_par);
    REPORT(s_CAApred);
    REPORT(s_CALpred);
    REPORT(s_vul);
    REPORT(s_L5);
    REPORT(s_LFS);
    REPORT(s_Vmaxlen);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif



