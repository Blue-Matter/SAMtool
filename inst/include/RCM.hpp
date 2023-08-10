
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

  DATA_IVECTOR(condition); // Indicates whether the model will condition on catch with estimated F (0), internally solve for F (1) or 'effort' (2)
  DATA_INTEGER(nll_C);    // Indicates whether there is a likelihood for the catch. TRUE when condition = 'catch' or 'effort' with nfleet > 1

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

  DATA_MATRIX(M_data);    // Natural mortality at age and year (n_y, n_age), overridden if there's a prior for log_M
  DATA_MATRIX(len_age);   // Length-at-age (n_y + 1, n_age)
  DATA_SCALAR(Linf);      // Linf
  DATA_SCALAR(K);
  DATA_SCALAR(t0);
  DATA_MATRIX(growth_time_f);
  DATA_MATRIX(growth_time_i);
  
  DATA_MATRIX(wt_len);
  
  DATA_MATRIX(SD_LAA);    // Length-at-age SD (n_y, n_age)
  DATA_MATRIX(wt);        // Weight-at-age stock (n_y + 1, n_age)
  DATA_MATRIX(mat);       // Maturity-at-age (n_y + 1, n_age) - only for matching survey selectivity to SSB
  DATA_MATRIX(fec);       // Fecundity-at-age (n_y + 1, n_age) - product of maturity and spawning output
  
  DATA_ARRAY(wt_c);       // Weight-at-age for fishery (n_y + 1, n_age, nfleet)

  DATA_IVECTOR(vul_type); // Integer vector indicating whether free (-2), logistic (-1), or dome vul (0) is used
  DATA_IVECTOR(ivul_type); // Same but for surveys, but can also mirror to B (-4), SSB (-3), or fleet (>0)
  DATA_IVECTOR(abs_I);    // Boolean, whether index is an absolute (fix q = 1) or relative terms (estimate q)
  DATA_IVECTOR(I_units);  // Boolean, whether index is biomass based (= 1) or abundance-based (0)

  DATA_MATRIX(age_error); // Ageing error matrix

  DATA_STRING(SR_type);   // String indicating whether Beverton-Holt ("BH"), Ricker, or "Mesnil-Rochet" stock-recruit is used
  DATA_MATRIX(LWT_fleet); // Likelihood weights for catch, C_eq, CAA, CAL, MS
  DATA_MATRIX(LWT_index); // Likelihood weights for the indices data
  DATA_STRING(comp_like); // Distribution for age/length comps (multinomial, lognormal, mvlogistic, dirmult1, dirmult2)

  DATA_SCALAR(max_F);     // Maximum F in the model
  DATA_SCALAR(rescale);   // R0 rescaler
  DATA_INTEGER(ageM);     // Age of maturity used for averaging E0 and EPR0 for stock-recruit relationship

  DATA_IVECTOR(yind_F);   // When condition(ff) = "catch", the year in F's are estimated and all other F parameters are deviations from this F
  DATA_INTEGER(n_itF);    // When condition(ff) = "catch2", the number of iterations for Newton-Raphson method to solve for F
  DATA_INTEGER(plusgroup) // Boolean, whether the maximum age in the plusgroup is modeled.
  
  DATA_IVECTOR(use_prior); // Boolean vector, whether to set a prior for R0, h, M, q (length of 3 + nsurvey)
  DATA_MATRIX(prior_dist); // Distribution of priors for R0, h, M, q (rows), columns indicate parameters of distribution calculated in R (see RCM_prior fn)
  
  DATA_INTEGER(nll_gr);    // Whether to ADREPORT annual likelihoods
  
  DATA_IMATRIX(est_vul);  // Indicates if fishery selectivity parameters are estimated. Then, will add vague priors to aid convergence
  DATA_IMATRIX(est_ivul); // Same but for index selectivity parameters
  
  DATA_IVECTOR(est_early_rec_dev); // Indicates years in which log_early_rec_dev are estimated. Then, lognormal bias correction estimates are added..
  DATA_IVECTOR(est_rec_dev); // Indicates years in which log_rec_dev are estimated.
  
  DATA_INTEGER(sim_process_error); // Whether to simulate process error (when using the TMB SIMULATE module) 
  DATA_SCALAR(spawn_time_frac);    // Fraction of year when spawning occurs for calculating spawning biomass
  DATA_IVECTOR(est_q);             // Whether to estimate index q (TRUE), otherwise solved analytically (FALSE)
  
  PARAMETER(R0x);                       // Unfished recruitment
  PARAMETER(transformed_h);             // Steepness
  PARAMETER_VECTOR(MR_SRR);             // Mesnil-Rochet parameters: (1) logit_Ehinge_E0, (2) log_delta
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
  
  PARAMETER_MATRIX(log_compf);          // CAA/CAL dispersion parameter
  PARAMETER_MATRIX(log_compi);          // IAA/IAL dispersion parameter
  
  PARAMETER_VECTOR(log_q);              // Index q

  int nlbin = lbinmid.size();
  Type R0 = 0;
  Type h = 0;
  if(SR_type == "BH") {
    h = 0.8 * invlogit(transformed_h) + 0.2;
    R0 = exp(R0x)/rescale;
  } else if(SR_type == "Ricker") {
    h = exp(transformed_h) + 0.2;
    R0 = exp(R0x)/rescale;
  }
  
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
  
  // Annual probability of length-at-age (PLA) - do check if needed
  vector<array<Type> > PLAf(n_y+1);
  vector<array<Type> > PLAi(n_y);
  
  int do_len = 0;
  if(CAL_n.sum() > 0 || IAL_n.sum() > 0 || (msize_type == "length" && !R_IsNA(asDouble(msize.sum())))) do_len += 1;
  for(int ff=0;ff<nfleet;ff++) if(vul_type(ff) == 0 || vul_type(ff) == -1) do_len += 1;
  for(int sur=0;sur<nsurvey;sur++) if(ivul_type(sur) == 0 || ivul_type(sur) == -1) do_len += 1;
  
  if (do_len > 0) {
    array<Type> Len_age_f = calculate_LAA(Linf, K, t0, growth_time_f, n_y+1, n_age, nfleet);
    array<Type> SD_LAA_f = calculate_SD_LAA(Len_age_f, Type(1.357), Type(0.001), n_y+1, n_age, nfleet);
    for(int y=0;y<=n_y;y++) {
      PLAf(y) = generate_PLA(lbin, Len_age_f, SD_LAA_f, n_age, nlbin, nfleet, y);
    }
    
    array<Type> Len_age_i = calculate_LAA(Linf, K, t0, growth_time_i, n_y, n_age, nsurvey);
    array<Type> SD_LAA_i = calculate_SD_LAA(Len_age_i, Type(1.357), Type(0.001), n_y, n_age, nsurvey);
    for(int y=0;y<n_y;y++) {
      PLAi(y) = generate_PLA(lbin, Len_age_i, SD_LAA_i, n_age, nlbin, nsurvey, y);
    }
  //  for(int y=0;y<=n_y;y++) PLA(y) = generate_PLA(lbin, len_age, SD_LAA, n_age, nlbin, y);
  }
  
  
  
  // Vulnerability (age or length-based) and F parameters
  vector<Type> LFS(nsel_block);
  vector<Type> L5(nsel_block);
  vector<Type> Vmaxlen(nsel_block);
  matrix<Type> vul_len(nlbin,nsel_block);
  vul_len.setZero();
  array<Type> vul = calc_vul(vul_par, vul_type, lbinmid, n_y+1, n_age, PLAf, LFS, L5, Vmaxlen, Linf, nfleet, sel_block, nsel_block, vul_len, prior, est_vul);

  vector<Type> q_effort(nfleet);
  vector<Type> CV_msize(nfleet);
  vector<Type> F_equilibrium(nfleet);
  F_equilibrium.setZero();

  matrix<Type> F(n_y,nfleet);
  matrix<Type> Z = M; // F to be added later during loop over years

  for(int ff=0;ff<nfleet;ff++) {
    CV_msize(ff) = exp(log_CV_msize(ff));
    q_effort(ff) = exp(log_q_effort(ff));
    if(condition(ff) != 2 && C_eq(ff)>0) F_equilibrium(ff) = exp(log_F_equilibrium(ff)); // catch/catch2
    if(condition(ff) == 2 && E_eq(ff)>0) F_equilibrium(ff) = q_effort(ff) * E_eq(ff);    // effort
    if(condition(ff) == 0) { //catch - set up F in midpoint of time series
      Type tmp = max_F - exp(log_F_dev(yind_F(ff),ff));
      F(yind_F(ff),ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), exp(log_F_dev(yind_F(ff),ff)));
    }
  }

  ////// Equilibrium reference points and per-recruit quantities - calculate annually
  vector<vector<Type> > NPR_unfished(n_y);
  vector<Type> EPR0(n_y);
  Type EPR0_SR = 0;
  for(int y=0;y<n_y;y++) {
    NPR_unfished(y) = calc_NPR0(M, n_age, y, plusgroup, spawn_time_frac);
    EPR0(y) = sum_EPR(NPR_unfished(y), fec, n_age, y);
    if(y <= ageM) EPR0_SR += EPR0(y);
  }
  EPR0_SR /= Type(ageM + 1); // Unfished replacement line for SRR is the average across the first ageM years
  Type E0_SR = 0;
  
  Type CR_SR = 0;
  Type Arec = 0;
  Type Brec = 0;
  Type MRRmax = 0;
  Type MRhinge = 0;
  Type MRgamma = 0;
  if(SR_type == "BH") {
    E0_SR = R0 * EPR0_SR;
    CR_SR = 4 *h;
    CR_SR /= 1-h;
    Brec = 5*h - 1;
    Brec /= (1-h);
    Brec /= E0_SR;
    Arec = CR_SR/EPR0_SR;
  } else if(SR_type == "Ricker") {
    E0_SR = R0 * EPR0_SR;
    CR_SR = pow(5*h, 1.25);
    Brec = 1.25;
    Brec *= log(5*h);
    Brec /= E0_SR;
    Arec = CR_SR/EPR0_SR;
  } else if(SR_type == "Mesnil-Rochet") {
    MRRmax = exp(R0x)/rescale;
    MRhinge = invlogit(MR_SRR(0)) * MRRmax * EPR0_SR;
    MRgamma = exp(MR_SRR(1));
    
    R0 = MesnilRochet_SR(EPR0_SR, MRgamma, MRRmax, MRhinge, 0);
    E0_SR = R0 * EPR0_SR;
    
    Type MR_K = pow(MRhinge * MRhinge + 0.25 * MRgamma * MRgamma, 0.5);
    Type MR_beta = R0/(MRhinge + MR_K);
    
    CR_SR = 2 * MR_beta * EPR0_SR;
  }
  
  ////// During time series year = 1, 2, ..., n_y
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
  vector<Type> Rec_dev(n_y);
  vector<Type> R_early(n_age-1);
  vector<Type> Rec_dev_early(n_age-1);
  matrix<Type> VB(n_y+1, nfleet);   // Vulnerable biomass at year
  vector<Type> B(n_y+1);            // Total biomass at year
  vector<Type> E(n_y+1);            // Spawning biomass at year
  
  vector<Type> log_rec_dev_sim = log_rec_dev;
  vector<Type> log_early_rec_dev_sim = log_early_rec_dev;
  
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
  
  Rec_dev.fill(1);
  Rec_dev_early.fill(1);
  
  log_rec_dev_sim.setZero();
  log_early_rec_dev_sim.setZero();
  
  // Equilibrium quantities (leading into first year of model)
  vector<Type> NPR_equilibrium = calc_NPR(F_equilibrium, vul, nfleet, M, n_age, 0, plusgroup, Type(0));
  vector<Type> NPR_equilibrium_spawn = calc_NPR(F_equilibrium, vul, nfleet, M, n_age, 0, plusgroup, spawn_time_frac);
  Type EPR_eq = sum_EPR(NPR_equilibrium_spawn, fec, n_age, 0);
  Type R_eq;

  if(SR_type == "BH") {
    R_eq = Arec * EPR_eq - 1;
    R_eq /= Brec * EPR_eq;
  } else if(SR_type == "Ricker") {
    R_eq = log(Arec * EPR_eq);
    R_eq /= Brec * EPR_eq;
  } else {
    R_eq = MesnilRochet_SR(EPR_eq, MRgamma, MRRmax, MRhinge, 0);
  }
  
  R(0) = R_eq;
  if(est_rec_dev(0)) {
    Rec_dev(0) = exp(log_rec_dev(0) - 0.5 * tau * tau);
    SIMULATE if(sim_process_error) {
      log_rec_dev_sim(0) = rnorm(log_rec_dev(0), tau);
      Rec_dev(0) = exp(log_rec_dev_sim(0) - 0.5 * tau * tau);
    }
    R(0) *= Rec_dev(0);
  }
  
  for(int a=0;a<n_age;a++) {
    if(a == 0) {
      N(0,a) = R(0) * NPR_equilibrium(a);
    } else {
      R_early(a-1) = R_eq;
      if(est_early_rec_dev(a-1)) {
        Rec_dev_early(a-1) = exp(log_early_rec_dev(a-1) - 0.5 * tau * tau);
        SIMULATE if(sim_process_error) {
          log_early_rec_dev_sim(a-1) = rnorm(log_early_rec_dev(a-1), tau);
          Rec_dev_early(a-1) = exp(log_early_rec_dev_sim(a-1) - 0.5 * tau * tau);
        }
        R_early(a-1) *= Rec_dev_early(a-1);
      }
      N(0,a) = R_early(a-1) * NPR_equilibrium(a);
    }
    
    Type Z_eq = M(0,a);
    for(int ff=0;ff<nfleet;ff++) Z_eq += vul(0,a,ff) * F_equilibrium(ff);
    for(int ff=0;ff<nfleet;ff++) {
      C_eq_pred(ff) += vul(0,a,ff) * F_equilibrium(ff) * wt_c(0,a,ff) * N(0,a) * (1 - exp(-Z_eq)) / Z_eq;
      VB(0,ff) += N(0,a) * wt_c(0,a,ff) * vul(0,a,ff);
    }
    
    B(0) += N(0,a) * wt(0,a);
  }

  // Loop over all other years
  for(int y=0;y<n_y;y++) {
    // Calculate this year's fleet F
    for(int ff=0;ff<nfleet;ff++) {
      if(condition(ff) == 0) { // catch
        if(y != yind_F(ff)) {
          Type tmp = max_F - F(yind_F(ff),ff) * exp(log_F_dev(y,ff)); // annual F as deviation from F in middle of time series
          F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), F(yind_F(ff),ff) * exp(log_F_dev(y,ff)));
        }
      } else if(condition(ff) == 1) { //catch2 - do not use due to wt_c
        F.row(y) = Newton_F(C_hist, N, M, wt, VB, vul, max_F, y, n_age, nfleet, n_itF, penalty);
      } else { //effort
        Type tmp = max_F - q_effort(ff) * E_hist(y,ff);
        F(y,ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), q_effort(ff) * E_hist(y,ff));
      }
    }
    
    // Calculate this year's total F, Z, SSB (ex. age-0)
    for(int a=0;a<n_age;a++) {
      for(int ff=0;ff<nfleet;ff++) Z(y,a) += vul(y,a,ff) * F(y,ff);
      if(a>0) E(y) += N(y,a) * exp(-spawn_time_frac * Z(y,a)) * fec(y,a);
    }
    
    // Calculate this year's recruitment and biomass
    if(y>0) {
      if(SR_type == "BH") {
        R(y) = BH_SR(E(y), h, R0, E0_SR);
      } else if(SR_type == "Ricker") {
        R(y) = Ricker_SR(E(y), h, R0, E0_SR);
      } else { // Mesnil-Rochet
        R(y) = MesnilRochet_SR(E(y), MRgamma, MRRmax, MRhinge);
      }
      
      if(est_rec_dev(y)) {
        Rec_dev(y) = exp(log_rec_dev(y) - 0.5 * tau * tau);
        SIMULATE if(sim_process_error) {
          log_rec_dev_sim(y) = rnorm(log_rec_dev(y), tau);
          Rec_dev(y) = exp(log_rec_dev_sim(y) - 0.5 * tau * tau);
        }
        R(y) *= Rec_dev(y);
      }
      N(y,0) = R(y);
      
      for(int a=0;a<n_age;a++) {
        B(y) += N(y,a) * wt(y,a);
        for(int ff=0;ff<nfleet;ff++) VB(y,ff) += vul(y,a,ff) * N(y,a) * wt_c(y,a,ff);
      }
    }
    
    // Calculate this year's CAA, catch, CAL, mean size, then next year's abundance
    for(int a=0;a<n_age;a++) {
      for(int ff=0;ff<nfleet;ff++) {
        CAAtrue(y,a,ff) = vul(y,a,ff) * F(y,ff) * N(y,a) * (1 - exp(-Z(y,a))) / Z(y,a);
        CN(y,ff) += CAAtrue(y,a,ff);
        Cpred(y,ff) += CAAtrue(y,a,ff) * wt_c(y,a,ff);
        for(int aa=0;aa<n_age;aa++) CAApred(y,aa,ff) += CAAtrue(y,a,ff) * age_error(a,aa); // a = true, aa = observed ages
        
        //if (CAL_hist.col(ff).sum() > 0) {
          for(int len=0;len<nlbin;len++) {
            CALpred(y,len,ff) += CAAtrue(y,a,ff) * PLAf(y)(a,len,ff);
            //Cpred(y,ff) += CALpred(y,len,ff) * wt_len(y,len);
          }
        //}
        if (msize_type == "length" && !R_IsNA(asDouble(msize.col(ff).sum())) && msize.col(ff).sum() > 0) {
          for(int len=0;len<nlbin;len++) MLpred(y,ff) += CAAtrue(y,a,ff) * PLAf(y)(a,len,ff) * lbinmid(len);
        }
      }
      
      if(a<n_age-1) N(y+1,a+1) = N(y,a) * exp(-Z(y,a));
      if(plusgroup && a==n_age-1) N(y+1,a) += N(y,a) * exp(-Z(y,a));
    }
    if (msize_type == "length") {
      for(int ff=0;ff<nfleet;ff++) {
        if (!R_IsNA(asDouble(msize.col(ff).sum())) && msize.col(ff).sum() > 0) {
          MLpred(y,ff) /= CN(y,ff);
        }
      }
    } else { //if(msize_type == "weight") 
      for(int ff=0;ff<nfleet;ff++) MWpred(y,ff) = Cpred(y,ff)/CN(y,ff);
    }
  }
  
  // Biomass at beginning of n_y + 1
  for(int a=1;a<n_age;a++) E(n_y) += N(n_y,a) * fec(n_y,a);
  
  if(spawn_time_frac > 0) { // Should work properly since spawn_time_frac is identified as DATA_SCALAR
    R(n_y) = R(n_y-1);
  } else if(SR_type == "BH") {
    R(n_y) = BH_SR(E(n_y), h, R0, E0_SR);
  } else if(SR_type == "Ricker") {
    R(n_y) = Ricker_SR(E(n_y), h, R0, E0_SR);
  } else { // Mesnil-Rochet
    R(n_y) = MesnilRochet_SR(E(n_y), MRgamma, MRRmax, MRhinge);
  }
  N(n_y,0) = R(n_y);
  for(int a=0;a<n_age;a++) {
    B(n_y) += N(n_y,a) * wt(n_y,a);
    for(int ff=0;ff<nfleet;ff++) VB(n_y,ff) += vul(n_y,a,ff) * N(n_y,a) * wt_c(n_y,a);
  }

  // Calculate for surveys: q, selectivity, and age/length comps
  vector<Type> iLFS(nsurvey);
  vector<Type> iL5(nsurvey);
  vector<Type> iVmaxlen(nsurvey);
  matrix<Type> ivul_len(nlbin, nsurvey);

  array<Type> IAAtrue(n_y, n_age, nsurvey); // True abundance at age vulnerable to survey
  array<Type> IAApred(n_y, n_age, nsurvey); // Predicted abundance (after ageing error) at age vulnerable to survey
  array<Type> IALpred(n_y, nlbin, nsurvey); // Abundance at length vulnerable to survey
  matrix<Type> IN(n_y, nsurvey);            // Total abundance vulnerable to the survey
  matrix<Type> Itot(n_y, nsurvey);          // Ipred/q - biomass or abundance vulnerable to the survey
  
  ivul_len.setZero();
  
  IAApred.setZero();
  IALpred.setZero();
  IN.setZero();
  Itot.setZero();

  array<Type> ivul = calc_ivul(ivul_par, ivul_type, lbinmid, n_y, n_age, PLAi, iLFS, iL5, iVmaxlen, Linf, 
                               mat, vul, ivul_len, prior, est_ivul);
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      for(int a=0;a<n_age;a++) {
        IAAtrue(y,a,sur) = ivul(y,a,sur) * N(y,a) * exp(-Z(y,a) * growth_time_i(y,sur));
        IN(y,sur) += IAAtrue(y,a,sur);

        for(int aa=0;aa<n_age;aa++) IAApred(y,aa,sur) += IAAtrue(y,a,sur) * age_error(a,aa);

        //if(I_units(sur)) Itot(y,sur) += IAAtrue(y,a,sur) * wt(y,a); // Biomass vulnerable to survey
        //if(IAL_n.col(sur).sum() > 0) { // Predict survey length comps if there are data
          for(int len=0;len<nlbin;len++) {
            IALpred(y,len,sur) += IAAtrue(y,a,sur) * PLAi(y)(a,len,sur);
            if(I_units(sur)) Itot(y,sur) += IALpred(y,len,sur) * wt_len(y,len);
          }
        //}
      }
    }
    if(!I_units(sur)) Itot.col(sur) = IN.col(sur); // Abundance vulnerable to survey
    if(est_q(sur)) {
      q(sur) = exp(log_q(sur));
      for(int y=0;y<n_y;y++) Ipred(y,sur) = q(sur) * Itot(y,sur);
    } else {
      q(sur) = calc_q(I_hist, Itot, sur, sur, Ipred, abs_I, n_y); // This function updates Ipred, uses '&'
    }
  }

  // Calc likelihood and parameter prior
  prior -= calc_prior(use_prior, prior_dist, R0x, h, SR_type == "BH", log_M, q, rescale);
  
  // MR penalty if MRhinge < min(E)
  if(SR_type == "Mesnil-Rochet") {
    penalty -= CppAD::CondExpGt(MRhinge, min(E), Type(0), dnorm(log(MRhinge), log(min(E)), Type(2), true));
  }
  
  array<Type> nll_fleet(n_y,nfleet,5);
  array<Type> nll_index(n_y,nsurvey,3);
  Type nll_log_rec_dev = 0;

  nll_fleet.setZero();
  nll_index.setZero();

  for(int sur=0;sur<nsurvey;sur++) {
    for(int y=0;y<n_y;y++) {
      if(LWT_index(sur,0) > 0 && !R_IsNA(asDouble(I_hist(y,sur)))) {
        nll_index(y,sur,0) -= LWT_index(sur,0) * dnorm_(log(I_hist(y,sur)), log(Ipred(y,sur)), sigma_I(y,sur), true);
        SIMULATE {
          I_hist(y,sur) = exp(rnorm(log(Ipred(y,sur)), sigma_I(y,sur)));
        }
      }
      
      if(LWT_index(sur,1) > 0 && !R_IsNA(asDouble(IAA_n(y,sur))) && IAA_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_index(y,sur,1) -= LWT_index(sur,1) * comp_multinom(IAA_hist, IAApred, IN, IAA_n, y, n_age, sur);
        } else if(comp_like == "lognormal") {
          nll_index(y,sur,1) -= LWT_index(sur,1) * comp_lognorm(IAA_hist, IAApred, IN, y, n_age, sur);
        } else if(comp_like == "dirmult1") {
          nll_index(y,sur,1) -= LWT_index(sur,1) * comp_dirmult1(IAA_hist, IAApred, IN, IAA_n, exp(log_compi(sur,0)), y, n_age, sur);
        } else if(comp_like == "dirmult2") {
          nll_index(y,sur,1) -= LWT_index(sur,1) * comp_dirmult2(IAA_hist, IAApred, IN, IAA_n, exp(log_compi(sur,0)), y, n_age, sur);
        }
      }
      
      if(LWT_index(sur,2) > 0 && !R_IsNA(asDouble(IAL_n(y,sur))) && IAL_n(y,sur) > 0) {
        if(comp_like == "multinomial") {
          nll_index(y,sur,2) -= LWT_index(sur,2) * comp_multinom(IAL_hist, IALpred, IN, IAL_n, y, nlbin, sur);
        } else if(comp_like == "lognormal") {
          nll_index(y,sur,2) -= LWT_index(sur,2) * comp_lognorm(IAL_hist, IALpred, IN, y, nlbin, sur);
        } else if(comp_like == "dirmult1") {
          nll_index(y,sur,2) -= LWT_index(sur,2) * comp_dirmult1(IAL_hist, IALpred, IN, IAL_n, exp(log_compi(sur,1)), y, nlbin, sur);
        } else if(comp_like == "dirmult2") {
          nll_index(y,sur,2) -= LWT_index(sur,2) * comp_dirmult2(IAL_hist, IALpred, IN, IAL_n, exp(log_compi(sur,1)), y, nlbin, sur);
        }
      }
    }
  }

  for(int ff=0;ff<nfleet;ff++) {
    if(LWT_fleet(ff,1) > 0 && C_eq(ff) > 0) {
      nll_fleet(0,ff,1) -= LWT_fleet(ff,1) * dnorm_(log(C_eq(ff)), log(C_eq_pred(ff)), sigma_Ceq(ff), true);
      SIMULATE {
        C_eq(ff) = exp(rnorm(log(C_eq_pred(ff)), sigma_Ceq(ff)));
      }
    }
    
    for(int y=0;y<n_y;y++) {
      int check1 = (condition(ff) != 2) && (C_hist(y,ff) > 0);
      int check2 = (condition(ff) == 2) && (E_hist(y,ff) > 0);
      if(check1 || check2) {
        
        if(nll_C && LWT_fleet(ff,0) > 0 && !R_IsNA(asDouble(C_hist(y,ff))) && C_hist(y,ff) > 0) {
          nll_fleet(y,ff,0) -= LWT_fleet(ff,0) * dnorm_(log(C_hist(y,ff)), log(Cpred(y,ff)), sigma_C(y,ff), true);
          
          SIMULATE {
            C_hist(y,ff) = exp(rnorm(log(Cpred(y,ff)), sigma_C(y,ff)));
          }
        } else if(!R_IsNA(asDouble(C_hist(y,ff))) && C_hist(y,ff) > 0) {
          SIMULATE {
            C_hist(y,ff) = Cpred(y,ff);
          }
        }
        
        if(LWT_fleet(ff,2) > 0 && !R_IsNA(asDouble(CAA_n(y,ff))) && CAA_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_fleet(y,ff,2) -= LWT_fleet(ff,2) * comp_multinom(CAA_hist, CAApred, CN, CAA_n, y, n_age, ff);
          } else if(comp_like == "lognormal") {
            nll_fleet(y,ff,2) -= LWT_fleet(ff,2) * comp_lognorm(CAA_hist, CAApred, CN, y, n_age, ff);
          } else if(comp_like == "dirmult1") {
            nll_fleet(y,ff,2) -= LWT_fleet(ff,2) * comp_dirmult1(CAA_hist, CAApred, CN, CAA_n, exp(log_compf(ff,0)), y, n_age, ff);
          } else if(comp_like == "dirmult2") {
            nll_fleet(y,ff,2) -= LWT_fleet(ff,2) * comp_dirmult2(CAA_hist, CAApred, CN, CAA_n, exp(log_compf(ff,0)), y, n_age, ff);
          }
        }
        
        if(LWT_fleet(ff,3) > 0 && !R_IsNA(asDouble(CAL_n(y,ff))) && CAL_n(y,ff) > 0) {
          if(comp_like == "multinomial") {
            nll_fleet(y,ff,3) -= LWT_fleet(ff,3) * comp_multinom(CAL_hist, CALpred, CN, CAL_n, y, nlbin, ff);
          } else if(comp_like == "lognormal") {
            nll_fleet(y,ff,3) -= LWT_fleet(ff,3) * comp_lognorm(CAL_hist, CALpred, CN, y, nlbin, ff);
          } else if(comp_like == "dirmult1") {
            nll_fleet(y,ff,3) -= LWT_fleet(ff,3) * comp_dirmult1(CAL_hist, CALpred, CN, CAL_n, exp(log_compf(ff,1)), y, nlbin, ff);
          } else if(comp_like == "dirmult2") {
            nll_fleet(y,ff,3) -= LWT_fleet(ff,3) * comp_dirmult2(CAL_hist, CALpred, CN, CAL_n, exp(log_compf(ff,1)), y, nlbin, ff);
          }
        }

        if(LWT_fleet(ff,4) > 0 && !R_IsNA(asDouble(msize(y,ff))) && msize(y,ff) > 0) {
          if(msize_type == "length") {
            nll_fleet(y,ff,4) -= LWT_fleet(ff,4) * dnorm_(msize(y,ff), MLpred(y,ff), CV_msize(ff) * msize(y,ff), true);
            SIMULATE {
              msize(y,ff) = rnorm(MLpred(y,ff), CV_msize(ff) * msize(y,ff));
            }
          } else {
            nll_fleet(y,ff,4) -= LWT_fleet(ff,4) * dnorm_(msize(y,ff), MWpred(y,ff), CV_msize(ff) * msize(y,ff), true);
            SIMULATE {
              msize(y,ff) = rnorm(MWpred(y,ff), CV_msize(ff) * msize(y,ff));
            }
          }
        }
      }
    }
  }

  for(int y=0;y<n_y;y++) {
    if(est_rec_dev(y) == 1) nll_log_rec_dev -= dnorm_(log_rec_dev(y), Type(0), tau, true);
  }
  for(int a=0;a<n_age-1;a++) {
    if(est_early_rec_dev(a) == 1) nll_log_rec_dev -= dnorm_(log_early_rec_dev(a), Type(0), tau, true);
  }
  
  if(comp_like == "mvlogistic") {
    for(int sur=0;sur<nsurvey;sur++) {
      if(LWT_index(sur,0) > 0 && IAA_n.col(sur).sum() > 0) {
        nll_index(0,sur,1) -= LWT_index(sur,1) * comp_mvlogistic(IAA_hist, IAApred, IN, n_y, n_age, sur);
      } 
      if(LWT_index(sur,2) > 0 && IAL_n.col(sur).sum() > 0) {
        nll_index(0,sur,2) -= LWT_index(sur,2) * comp_mvlogistic(IAL_hist, IALpred, IN, n_y, nlbin, sur);
      }
    }
    
    for(int ff=0;ff<nfleet;ff++) {
      if(LWT_fleet(ff,2) > 0 && CAA_n.col(ff).sum() > 0) {
        nll_fleet(0,ff,2) -= LWT_fleet(ff,2) * comp_mvlogistic(CAA_hist, CAApred, CN, n_y, n_age, ff);
      }
      if(LWT_fleet(ff,3) > 0 && CAL_n.col(ff).sum() > 0) {
        nll_fleet(0,ff,3) -= LWT_fleet(ff,3) * comp_mvlogistic(CAL_hist, CALpred, CN, n_y, nlbin, ff);
      }
    }
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
  //if(condition == "effort") ADREPORT(q_effort);
  if(nll_index.col(0).sum() != 0) ADREPORT(q);

  REPORT(R0x);
  REPORT(transformed_h);
  
  if(SR_type == "Mesnil-Rochet") {
    REPORT(MR_SRR);
    REPORT(MRRmax);
    REPORT(MRhinge);
    REPORT(MRgamma);
    ADREPORT(MRRmax);
    ADREPORT(MRhinge);
    ADREPORT(MRgamma);
  }
  
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
  //if(condition == "catch") REPORT(log_F_dev);
  REPORT(F_equilibrium);
  REPORT(vul);
  if(vul_len.sum() > 0) REPORT(vul_len);
  REPORT(F);
  REPORT(Z);

  REPORT(NPR_unfished);
  REPORT(EPR0);
  
  if(SR_type != "Mesnil-Rochet") {
    REPORT(Arec);
    REPORT(Brec);
  }
  REPORT(E0_SR);
  REPORT(EPR0_SR);
  REPORT(CR_SR);

  REPORT(N);
  REPORT(CAApred);
  if(CALpred.sum() > 0) REPORT(CALpred);
  if(MLpred.sum() > 0) REPORT(MLpred);
  if(msize_type == "weight" && MWpred.sum() > 0) REPORT(MWpred);
  REPORT(CN);
  REPORT(Cpred);
  REPORT(Ipred);
  REPORT(R);
  REPORT(R_early);
  REPORT(VB);
  REPORT(B);
  REPORT(E);
  REPORT(Rec_dev);

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
  
  if(comp_like == "dirmult1" || comp_like == "dirmult2") {
    REPORT(log_compf);
    REPORT(log_compi);
  }

  if(age_error.trace() != Type(n_age)) {
    REPORT(CAAtrue);
    REPORT(IAAtrue);
  }

  if(nll_index.sum() != 0) {
    REPORT(ivul_par);
    REPORT(ivul);
    if(ivul_len.sum() > 0) REPORT(ivul_len);
    REPORT(iL5);
    REPORT(iLFS);
    REPORT(iVmaxlen);
    REPORT(IAApred);
  }
  if(IAL_n.sum() > 0) REPORT(IALpred);
  
  if(nll_gr) {
    ADREPORT(nll_fleet);
    ADREPORT(nll_index);
  }
  
  SIMULATE {
    REPORT(C_eq);
    REPORT(C_hist);
    REPORT(I_hist);
    REPORT(msize);
    
    REPORT(log_rec_dev_sim);
    REPORT(log_early_rec_dev_sim);
  }

  return nll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif



