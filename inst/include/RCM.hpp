
#ifndef RCM_hpp
#define RCM_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type RCM(objective_function<Type> *obj) {

  using namespace ns_RCM;

  DATA_MATRIX(C_hist);    // Total catch by year and fleet
  DATA_VECTOR(C_eq);      // Equilibrium catch by fleet
  DATA_MATRIX(sigma_C);   // Standard deviation of catch
  DATA_VECTOR(sigma_Ceq); // Standard deviation of equilibrium catch

  DATA_MATRIX(E_hist);    // Effort by year and fleet
  DATA_VECTOR(E_eq);      // Equilibrium effort by fleet

  DATA_STRING(condition); // Indicates whether the model will condition on 'effort', 'catch' (estimated F's), 'catch2' (internally optimize for F)
  DATA_INTEGER(nll_C);    // Indicates whether there is a likelihood for the catch. TRUE when (a) condition = 'catch' or (b) 'effort' with nfleet > 1

  DATA_MATRIX(I_hist);    // Index by year and survey
  DATA_MATRIX(sigma_I);   // Standard deviation of index by year and survey

  DATA_ARRAY(CAA_hist);   // Catch-at-age proportions by year, age, fleet
  DATA_MATRIX(CAA_n);     // Annual samples in CAA by year and fleet

  DATA_ARRAY(CAL_hist);   // Catch-at-length proportions by year, length_bin, fleet
  DATA_MATRIX(CAL_n);     // Annual samples in CAL by year and fleet

  DATA_ARRAY(IAA_hist);   // Index-at-age proportions by year, age, survey
  DATA_MATRIX(IAA_n);     // Annual samples in IAA by year and survey

  DATA_ARRAY(IAL_hist);   // Index-at-length proportions by year, length_bin, survey
  DATA_MATRIX(IAL_n);     // Annual samples in IAL by year and survey

  DATA_VECTOR(lbin);      // Vector of length bin boundaries
  DATA_VECTOR(lbinmid);   // Vector of length bin midpoints
  DATA_MATRIX(msize);      // Vector of annual mean size by year and fleet
  DATA_STRING(msize_type); // Whether the mean size is length or weight

  DATA_IMATRIX(sel_block); // Assigns selectivity by year and fleet
  DATA_INTEGER(nsel_block); // The number of selectivity "blocks"

  DATA_INTEGER(n_y);      // Number of years in model
  DATA_INTEGER(n_age);    // Number of age classes from 0 - maxage (plus-group)
  DATA_INTEGER(nfleet);   // Number of fleets
  DATA_INTEGER(nsurvey);  // Number of surveys

  DATA_MATRIX(M_data);    // Natural mortality at age and year, overridden if there's a prior for log_M
  DATA_MATRIX(len_age);   // Length-at-age
  DATA_SCALAR(Linf);      // Linf
  DATA_MATRIX(SD_LAA);    // Length-at-age SD (year x age)
  DATA_MATRIX(wt);        // Weight-at-age
  DATA_MATRIX(mat);       // Maturity-at-age at the beginning of the year

  DATA_IVECTOR(vul_type); // Integer vector indicating whether free (-2), logistic (-1), or dome vul (0) is used
  DATA_IVECTOR(ivul_type); // Same but for surveys, but can also mirror to B (-4), SSB (-3), or fleet (>0)
  DATA_IVECTOR(abs_I);    // Boolean, whether index is an absolute (fix q = 1) or relative terms (estimate q)
  DATA_IVECTOR(I_units);  // Boolean, whether index is biomass based (= 1) or abundance-based (0)

  DATA_MATRIX(age_error); // Ageing error matrix

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt or Ricker stock-recruit is used
  DATA_MATRIX(LWT_fleet); // LIkelihood weights for catch, C_eq, CAA, CAL, MS
  DATA_MATRIX(LWT_index); // Likelihood weights for the indices data
  DATA_STRING(comp_like); // Whether to use "multinomial" or "lognormal" distribution for age/lengthc comps

  DATA_SCALAR(max_F);     // Maximum F in the model
  DATA_SCALAR(rescale);   // R0 rescaler
  DATA_INTEGER(ageM);     // Age of maturity used for averaging E0 and EPR0

  DATA_IVECTOR(yind_F);   // When condition = "catch", the year in F's are estimated and all other F parameters are deviations from this F
  DATA_INTEGER(n_itF);    // When condition = "catch2", the number of iterations for Newton-Raphson method to solve for F
  DATA_INTEGER(plusgroup) // Boolean, whether the maximum age in the plusgroup is modeled.
  
  DATA_IVECTOR(use_prior); // Boolean vector, whether to set a prior for R0, h, M, q (length of 3 + nsurvey)
  DATA_MATRIX(prior_dist); // Distribution of priors for R0, h, M, q (rows), columns indicate parameters of distribution calculated in R (see RCM_prior fn)
  
  DATA_INTEGER(nll_gr);    // Whether to ADREPORT annual likelihoods
  DATA_IVECTOR(est_early_rec_dev); // Indicates years in which log_early_rec_dev are estimated. Then, lognormal bias correction estimates are added..
  DATA_IVECTOR(est_rec_dev); // Indicates years in which log_rec_dev are estimated.

  PARAMETER(R0x);                       // Unfished recruitment
  PARAMETER(transformed_h);             // Steepness
  PARAMETER(log_M);                     // Age and time constant M (only if there's a prior, then it will override M_data)
  PARAMETER_MATRIX(vul_par);            // Matrix of vul_par 3 rows and nsel_block columns
  PARAMETER_MATRIX(ivul_par);           // Matrix of index selectivity parameters, 3 rows and nsurvey columns
  PARAMETER_VECTOR(log_q_effort);       // log_q for F when condition = "effort"
  PARAMETER_MATRIX(log_F_dev);          // log_F_deviations when condition = "catch"
  PARAMETER_VECTOR(log_F_equilibrium);  // Equilibrium F by fleet when condition != "effort"

  PARAMETER_VECTOR(log_CV_msize);       // CV of mean size
  PARAMETER(log_tau);                   // Std. dev. of rec devs
  PARAMETER_VECTOR(log_early_rec_dev);  // Rec devs for first year abundance
  PARAMETER_VECTOR(log_rec_dev);        // Rec devs for all other years in the model

  int nlbin = lbinmid.size();

  Type R0 = exp(R0x)/rescale;
  Type h;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h);
  } else {
    h = exp(transformed_h);
  }
  h += 0.2;
  Type Mest = exp(log_M);
  matrix<Type> M(n_y, n_age);
  if(use_prior(2)) {
    M.fill(Mest);
  } else {
    M = M_data;
  }
  Type tau = exp(log_tau);
  
  Type penalty = 0;
  Type prior = 0;
  
  // Vulnerability (length-based) and F parameters
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
      Type tmp = max_F - exp(log_F_dev(yind_F(ff),ff));
      F(yind_F(ff),ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), exp(log_F_dev(yind_F(ff),ff)));
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
    NPR_unfished(y) = calc_NPR0(M, n_age, y, plusgroup);

    EPR0(y) = sum_EPR(NPR_unfished(y), wt, mat, n_age, y);
    E0(y) = R0 * EPR0(y);
    B0(y) = R0 * sum_BPR(NPR_unfished(y), wt, n_age, y);
    N0(y) = R0 * NPR_unfished(y).sum();

    if(y <= ageM) E0_SR += E0(y);
  }
  E0_SR /= Type(ageM + 1);
  Type EPR0_SR = E0_SR/R0;

  Type CR_SR, Brec;
  if(SR_type == "BH") {
    CR_SR = 4 *h;
    CR_SR /= 1-h;
    Brec = 5*h - 1;
    Brec /= (1-h);
  } else {
    CR_SR = pow(5*h, 1.25);
    Brec = 1.25;
    Brec *= log(5*h);
  }
  Type Arec = CR_SR/EPR0_SR;
  Brec /= E0_SR;
  
  ////// During time series year = 1, 2, ..., n_y
  vector<matrix<Type> > ALK(n_y);
  matrix<Type> N(n_y+1, n_age);

  vector<Type> C_eq_pred(nfleet);
  array<Type> CAAtrue(n_y, n_age, nfleet);   // Catch (in numbers) at year and age at the mid-point of the season
  array<Type> CAApred(n_y, n_age, nfleet);   // Catch (in numbers) at year and age at the mid-point of the season
  array<Type> CALpred(n_y, nlbin, nfleet);
  matrix<Type> MLpred(n_y, nfleet);
  matrix<Type> MWpred(n_y, nfleet);
  matrix<Type> CN(n_y, nfleet);             // Catch in numbers

  matrix<Type> Cpred(n_y, nfleet);
  matrix<Type> Ipred(n_y, nsurvey);          // Predicted index at year

  vector<Type> R(n_y+1);            // Recruitment at year
  vector<Type> R_early(n_age-1);
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
  vector<Type> NPR_equilibrium = calc_NPR(F_equilibrium, vul, nfleet, M, n_age, 0, plusgroup);
  Type EPR_eq = sum_EPR(NPR_equilibrium, wt, mat, n_age, 0);
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
    
    B(0) += N(0,a) * wt(0,a);
    E(0) += N(0,a) * wt(0,a) * mat(0,a);

    Type Z_eq = M(0,a);
    for(int ff=0;ff<nfleet;ff++) Z_eq += vul(0,a,ff) * F_equilibrium(ff);
    for(int ff=0;ff<nfleet;ff++) {
      C_eq_pred(ff) += vul(0,a,ff) * F_equilibrium(ff) * wt(0,a) * N(0,a) * (1 - exp(-Z_eq)) / Z_eq;
      VB(0,ff) += N(0,a) * wt(0,a) * vul(0,a,ff);
    }
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    // Calculate this year's age-length key (ALK) and F
    if(Type(n_age) != Linf) ALK(y) = generate_ALK(lbin, len_age, SD_LAA, n_age, nlbin, y);
    if(condition == "catch") {
      for(int ff=0;ff<nfleet;ff++) {
        if(y != yind_F(ff)) {
          Type tmp = max_F - F(yind_F(ff),ff) * exp(log_F_dev(y,ff));
          F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), F(yind_F(ff),ff) * exp(log_F_dev(y,ff)));
        }
      }
    } else if(condition == "catch2") {
      F.row(y) = Newton_F(C_hist, N, M, wt, VB, vul, max_F, y, n_age, nfleet, n_itF, penalty);
    } else {
      for(int ff=0;ff<nfleet;ff++) {
        Type tmp = max_F - q_effort(ff) * E_hist(y,ff);
        F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), q_effort(ff) * E_hist(y,ff));
      }
    }
    
    // Calculate this year's catch, CAA, CAL, mean size; and next year's abundance and SSB (ex. age-0)
    for(int a=0;a<n_age;a++) {
      for(int ff=0;ff<nfleet;ff++) Z(y,a) += vul(y,a,ff) * F(y,ff);
      if(a<n_age-1) N(y+1,a+1) = N(y,a) * exp(-Z(y,a));
      if(plusgroup && a==n_age-1) N(y+1,a) += N(y,a) * exp(-Z(y,a));

      for(int ff=0;ff<nfleet;ff++) {
        CAAtrue(y,a,ff) = vul(y,a,ff) * F(y,ff) * N(y,a) * (1 - exp(-Z(y,a))) / Z(y,a);
        CN(y,ff) += CAAtrue(y,a,ff);
        Cpred(y,ff) += CAAtrue(y,a,ff) * wt(y,a);

        if(Type(n_age) != Linf) {
          for(int len=0;len<nlbin;len++) {
            CALpred(y,len,ff) += CAAtrue(y,a,ff) * ALK(y)(a,len);
            MLpred(y,ff) += CAAtrue(y,a,ff) * ALK(y)(a,len) * lbinmid(len);
          }
        }
        for(int aa=0;aa<n_age;aa++) CAApred(y,aa,ff) += CAAtrue(y,a,ff) * age_error(a,aa); // a = true, aa = observed ages
      }
      if(a>0) E(y+1) += N(y+1,a) * wt(y+1,a) * mat(y+1,a);
    }
    if(Type(n_age) != Linf) for(int ff=0;ff<nfleet;ff++) MLpred(y,ff) /= CN(y,ff);
    if(msize_type == "weight") for(int ff=0;ff<nfleet;ff++) MWpred(y,ff) = Cpred(y,ff)/CN(y,ff);
    
    // Calc next year's recruitment, total biomass, and vulnerable biomass
    if(SR_type == "BH") {
      R(y+1) = BH_SR(E(y+1), h, R0, E0_SR);
    } else {
      R(y+1) = Ricker_SR(E(y+1), h, R0, E0_SR);
    }
    
    if(y<n_y-1 && est_rec_dev(y+1)) R(y+1) *= exp(log_rec_dev(y+1) - 0.5 * tau * tau);
    N(y+1,0) = R(y+1);
    
    for(int a=0;a<n_age;a++) {
      B(y+1) += N(y+1,a) * wt(y+1,a);
      for(int ff=0;ff<nfleet;ff++) VB(y+1,ff) += vul(y+1,a,ff) * N(y+1,a) * wt(y+1,a);
    }
  }

  // Calculate survey q, selectivity, and age/length comps
  vector<Type> iLFS(nsurvey);
  vector<Type> iL5(nsurvey);
  vector<Type> iVmaxlen(nsurvey);

  array<Type> IAAtrue(n_y, n_age, nsurvey); // True abundance at age vulnerable to survey
  array<Type> IAApred(n_y, n_age, nsurvey); // Predicted abundance (after ageing error) at age vulnerable to survey
  array<Type> IALpred(n_y, nlbin, nsurvey); // Abundance at length vulnerable to survey
  matrix<Type> IN(n_y, nsurvey);            // Total abundance vulnerable to the survey
  matrix<Type> Itot(n_y, nsurvey);          // Ipred/q - biomass or abundance vulnerable to the survey

  IAApred.setZero();
  IALpred.setZero();
  IN.setZero();
  Itot.setZero();

  array<Type> ivul = calc_ivul(ivul_par, ivul_type, len_age, iLFS, iL5, iVmaxlen, Linf, mat, vul, prior);
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      for(int a=0;a<n_age;a++) {
        IAAtrue(y,a,sur) = ivul(y,a,sur) * N(y,a);
        IN(y,sur) += IAAtrue(y,a,sur);

        for(int aa=0;aa<n_age;aa++) IAApred(y,aa,sur) += IAAtrue(y,a,sur) * age_error(a,aa);

        if(I_units(sur)) Itot(y,sur) += IAAtrue(y,a,sur) * wt(y,a); // Biomass vulnerable to survey
        if(Type(n_age) != Linf && IAL_n.col(sur).sum() > 0) { // Predict survey length comps if there are data
          for(int len=0;len<nlbin;len++) IALpred(y,len,sur) += IAAtrue(y,a,sur) * ALK(y)(a,len);
        }
      }
    }
    if(!I_units(sur)) Itot.col(sur) = IN.col(sur); // Abundance vulnerable to survey
    q(sur) = calc_q(I_hist, Itot, sur, sur, Ipred, abs_I, n_y); // This function updates Ipred
  }

  // Calc likelihood and parameter prior
  prior -= calc_prior(use_prior, prior_dist, R0, h, SR_type == "BH", log_M, q);
  
  array<Type> nll_fleet(n_y,nfleet,5);
  array<Type> nll_index(n_y,nsurvey,3);
  Type nll_log_rec_dev = 0;

  nll_fleet.setZero();
  nll_index.setZero();

  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      if(LWT_index(sur,0) > 0 && !R_IsNA(asDouble(I_hist(y,sur)))) {
        nll_index(y,sur,0) -= LWT_index(sur,0) * dnorm_(log(I_hist(y,sur)), log(Ipred(y,sur)), sigma_I(y,sur), true);
      }
      
      if(LWT_index(sur,1) > 0 && !R_IsNA(asDouble(IAA_n(y,sur))) && IAA_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_index(y,sur,1) -= LWT_index(sur,1) * comp_multinom(IAA_hist, IAApred, IN, IAA_n, y, n_age, sur);
        } else {
          nll_index(y,sur,1) -= LWT_index(sur,1) * comp_lognorm(IAA_hist, IAApred, IN, y, n_age, sur);
        }
      }
      
      if(LWT_index(sur,2) > 0 && !R_IsNA(asDouble(IAL_n(y,sur))) && IAL_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_index(y,sur,2) -= LWT_index(sur,2) * comp_multinom(IAL_hist, IALpred, IN, IAL_n, y, nlbin, sur);
        } else {
          nll_index(y,sur,2) -= LWT_index(sur,2) * comp_lognorm(IAL_hist, IALpred, IN, y, nlbin, sur);
        }
      }
    }
  }

  for(int ff=0;ff<nfleet;ff++) {
    if(LWT_fleet(ff,1) > 0 && C_eq(ff) > 0) {
      nll_fleet(0,ff,1) -= LWT_fleet(ff,1) * dnorm_(log(C_eq(ff)), log(C_eq_pred(ff)), sigma_Ceq(ff), true);
    }
    
    for(int y=0;y<n_y;y++) {
      int check1 = (condition != "effort") && (C_hist(y,ff) > 0);
      int check2 = (condition == "effort") && (E_hist(y,ff) > 0);
      if(check1 || check2) {
        
        if(nll_C && LWT_fleet(ff,0) > 0 && !R_IsNA(asDouble(C_hist(y,ff))) && C_hist(y,ff) > 0) {
          nll_fleet(y,ff,0) -= LWT_fleet(ff,0) * dnorm_(log(C_hist(y,ff)), log(Cpred(y,ff)), sigma_C(y,ff), true);
        }
        
        if(LWT_fleet(ff,2) > 0 && !R_IsNA(asDouble(CAA_n(y,ff))) && CAA_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_fleet(y,ff,2) -= LWT_fleet(ff,2) * comp_multinom(CAA_hist, CAApred, CN, CAA_n, y, n_age, ff);
          } else {
            nll_fleet(y,ff,2) -= LWT_fleet(ff,2) * comp_lognorm(CAA_hist, CAApred, CN, y, n_age, ff);
          }
        }
        
        if(LWT_fleet(ff,3) > 0 && !R_IsNA(asDouble(CAL_n(y,ff))) && CAL_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_fleet(y,ff,3) -= LWT_fleet(ff,3) * comp_multinom(CAL_hist, CALpred, CN, CAL_n, y, nlbin, ff);
          } else {
            nll_fleet(y,ff,3) -= LWT_fleet(ff,3) * comp_lognorm(CAL_hist, CALpred, CN, y, nlbin, ff);
          }
        }

        if(LWT_fleet(ff,4) > 0 && !R_IsNA(asDouble(msize(y,ff))) && msize(y,ff) > 0) {
          if(msize_type == "length") {
            nll_fleet(y,ff,4) -= LWT_fleet(ff,4) * dnorm_(msize(y,ff), MLpred(y,ff), CV_msize(ff) * msize(y,ff), true);
          } else {
            nll_fleet(y,ff,4) -= LWT_fleet(ff,4) * dnorm_(msize(y,ff), MWpred(y,ff), CV_msize(ff) * msize(y,ff), true);
          }
        }
      }
    }
  }

  for(int y=0;y<n_y;y++) {
    if(est_rec_dev(y)) nll_log_rec_dev -= dnorm_(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<n_age-1;a++) {
    if(est_early_rec_dev(a)) nll_log_rec_dev -= dnorm_(log_early_rec_dev(a), Type(0), tau, true);
  }

  Type nll = nll_fleet.sum() + nll_index.sum();
  nll += nll_log_rec_dev + penalty + prior;

  if(CppAD::Variable(R0x)) ADREPORT(R0);
  ADREPORT(h);
  if(use_prior(2)) {
    ADREPORT(Mest);
    REPORT(Mest);
  }
  if(CppAD::Variable(log_tau)) ADREPORT(tau);
  if(condition == "effort") ADREPORT(q_effort);
  if(nll_index.col(0).sum() != 0) ADREPORT(q);

  REPORT(R0x);
  REPORT(transformed_h);
  REPORT(vul_par);
  REPORT(LFS);
  REPORT(L5);
  REPORT(Vmaxlen);
  REPORT(log_q_effort);
  REPORT(log_F_equilibrium);

  REPORT(log_early_rec_dev);
  REPORT(log_rec_dev);

  REPORT(R0);
  REPORT(h);
  REPORT(tau);
  if(nll_fleet.col(4).sum() != 0) REPORT(CV_msize);
  if(condition == "catch") REPORT(log_F_dev);
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
  REPORT(CR_SR);

  if(nll_fleet.col(3).sum() != 0 || nll_index.col(2).sum() != 0 || ((nll_fleet.col(4).sum() != 0) & (msize_type == "length"))) REPORT(ALK);
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
  if(nll_index.col(0).sum() != 0) REPORT(q);

  REPORT(nll_fleet);
  REPORT(nll_index);
  REPORT(nll_log_rec_dev);

  REPORT(nll);
  REPORT(penalty);
  REPORT(prior);

  if(age_error.trace() != Type(n_age)) {
    REPORT(CAAtrue);
    REPORT(IAAtrue);
  }

  if(nll_index.sum() != 0) {
    REPORT(ivul_par);
    REPORT(IAApred);
    REPORT(IALpred);
    REPORT(ivul);
    REPORT(iL5);
    REPORT(iLFS);
    REPORT(iVmaxlen);
  }
  
  if(nll_gr) {
    ADREPORT(nll_fleet);
    ADREPORT(nll_index);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif



