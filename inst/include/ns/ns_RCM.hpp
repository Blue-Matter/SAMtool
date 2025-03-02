
namespace ns_RCM {

using namespace ns_SCA;

template <class Type>
Type log2(Type x) {
  return log(x)/log(Type(2));
}

template<class Type>
matrix<Type> generate_PLA(vector<Type> lbin, matrix<Type> len_age, matrix<Type> SD_LAA, 
                          int n_age, int nlbin, int y) {
  matrix<Type> PLA(n_age, nlbin);
  for(int a=0;a<n_age;a++) {
    for(int j=0;j<nlbin;j++) {
      if(j==nlbin-1) {
        PLA(a,j) = 1 - pnorm(lbin(j), len_age(y,a), SD_LAA(y,a));
      } else {
        PLA(a,j) = pnorm(lbin(j+1), len_age(y,a), SD_LAA(y,a));
        if(j>0) PLA(a,j) -= pnorm(lbin(j), len_age(y,a), SD_LAA(y,a));
      }
    }
  }
  return PLA;
}

template<class Type>
vector<Type> calc_NPR0(matrix<Type> M, int n_age, int y, int plusgroup, Type spawn_time_frac = 0) {
  vector<Type> NPR(n_age);
  vector<Type> out(n_age);
  NPR(0) = 1;
  for(int a=1;a<n_age;a++) NPR(a) = NPR(a-1) * exp(-M(y,a-1));
  if(plusgroup) NPR(n_age-1) /= 1 - exp(-M(y,n_age-1));
  for(int a=0;a<n_age;a++) out(a) = NPR(a) * exp(-spawn_time_frac * M(y,a));
  return out;
}


template<class Type>
vector<Type> calc_NPR(vector<Type> F, array<Type> vul, int nfleet, matrix<Type> M, int n_age, int y, int plusgroup,
                      Type spawn_time_frac = 0) {
  vector<Type> NPR(n_age);
  vector<Type> Z = M.row(y);
  vector<Type> out(n_age);
  NPR(0) = 1;
  for(int a=0;a<n_age;a++) {
    for(int ff=0;ff<nfleet;ff++) Z(a) += vul(y,a,ff) * F(ff);
    if(a > 0) NPR(a) = NPR(a-1) * exp(-Z(a-1));
  }
  if(plusgroup) NPR(n_age-1) /= 1 - exp(-Z(n_age-1));
  for(int a=0;a<n_age;a++) out(a) = NPR(a) * exp(-spawn_time_frac * Z(a));
  return out;
}

template<class Type>
Type sum_EPR(vector<Type> NPR, matrix<Type> fec, int n_age, int y) {
  Type EPR = 0.;
  for(int a=0;a<n_age;a++) EPR += NPR(a) * fec(y,a);
  return EPR;
}


template<class Type>
array<Type> calc_vul(matrix<Type> vul_par, vector<int> vul_type, vector<Type> lbinmid, int n_y, int n_age,
                     vector<matrix<Type> > PLA, vector<Type> &LFS, vector<Type> &L5,
                     vector<Type> &Vmaxlen, Type Linf, int nfleet, matrix<int> sel_block, int nsel_block, 
                     matrix<Type> &vul_len, Type &prior, matrix<int> est_vul) {
  
  int nlbin = lbinmid.size();
  array<Type> vul(n_y, n_age, nfleet); // Corresponding age based selectivity
  vul.setZero();
  vector<Type> sls(nsel_block);
  vector<Type> srs(nsel_block);
  Type maxage = Type(n_age) - 1;
  
  for(int b=0;b<nsel_block;b++) { // Parameters for sel_block - don't do anything for age-specific index
    if(vul_type(b) == -2) { // Free parameters - adding priors only
      for(int a=0;a<n_age;a++) {
        if(est_vul(a,b)) {
          Type v = invlogit(vul_par(a,b));
          prior -= dbeta_(v, Type(1.01), Type(1.01), true) + log(v - v * v);
        }
      }
    } else {
      // Dome or logistic sel
      if(est_vul(1,b)) prior -= dnorm_(vul_par(1,b), Type(0), Type(3), true);
      
      Type ilogit_x = invlogit(vul_par(0,b));
      if(est_vul(0,b)) {
        Type jac = ilogit_x - ilogit_x * ilogit_x;
        prior -= dbeta_(ilogit_x, Type(1.01), Type(1.01), true) + log(jac);
        //prior -= dnorm_(vul_par(0,b), Type(0), Type(3), true);
      }
      
      if(vul_type(b) == 0 || vul_type(b) == - 1) {
        LFS(b) = ilogit_x * Linf;
      } else {
        LFS(b) = ilogit_x * maxage;
      }
      L5(b) = LFS(b) - exp(vul_par(1,b));
      sls(b) = (LFS(b) - L5(b))/pow(-log2(0.05), 0.5);
      
      if(vul_type(b) == -1 || vul_type(b) == -6) { // Logistic
        Vmaxlen(b) = 1;
      } else { // Dome
        Vmaxlen(b) = invlogit(vul_par(2,b));
        
        if(vul_type(b) == 0) {
          srs(b) = (Linf - LFS(b))/pow(-log2(Vmaxlen(b)), 0.5);
        } else {
          srs(b) = (maxage - LFS(b))/pow(-log2(Vmaxlen(b)), 0.5);
        }
        if(est_vul(2,b)) {
          Type jac = Vmaxlen(b) - Vmaxlen(b) * Vmaxlen(b);
          prior -= dbeta_(Vmaxlen(b), Type(1.01), Type(1.01), true) + log(jac);
        }
      }
      
      if(vul_type(b) == 0 || vul_type(b) == -1) { // Calculate length-based sel
        for(int j=0;j<nlbin;j++) {
          Type lo = pow(2, -((lbinmid(j) - LFS(b))/sls(b) * (lbinmid(j) - LFS(b))/sls(b)));
          Type hi;
          if(vul_type(b) == -1) {
            hi = 1;
          } else {
            hi = pow(2, -((lbinmid(j) - LFS(b))/srs(b) * (lbinmid(j) - LFS(b))/srs(b)));
          }
          vul_len(j,b) = CppAD::CondExpLt(lbinmid(j), LFS(b), lo, hi);
        }
      }
    }
  }
  
  for(int ff=0;ff<nfleet;ff++) { // Assign to fleet
    for(int y=0;y<n_y;y++) {
      int vul_ind = sel_block(y,ff) - 1;
      if(vul_type(vul_ind) == 0 || vul_type(vul_ind) == -1) { // Convert length sel to age sel
        vector<Type> vul_age(n_age); // Ensure max(vul) = 1
        vul_age.setZero();
        for(int a=0;a<n_age;a++) {
          for(int j=0;j<nlbin;j++) vul_age(a) += PLA(y)(a,j) * vul_len(j,vul_ind);
        }
        
        Type max_vul = max(vul_age);
        for(int a=0;a<n_age;a++) vul(y,a,ff) = vul_age(a)/max_vul;
      } else if(vul_type(vul_ind) == -5 || vul_type(vul_ind) == -6) {  // Dome or logistic, age
        vector<Type> vul_age(n_age); // Ensure max(vul) = 1
        
        for(int a=0;a<n_age;a++) { // Calculate age-based sel
          Type aa = Type(a);
          Type lo = pow(2, -((aa - LFS(vul_ind))/sls(vul_ind) * (aa - LFS(vul_ind))/sls(vul_ind)));
          Type hi;
          if(vul_type(vul_ind) == -6) {
            hi = 1;
          } else {
            hi = pow(2, -((aa - LFS(vul_ind))/srs(vul_ind) * (aa - LFS(vul_ind))/srs(vul_ind)));
          }
          vul_age(a) = CppAD::CondExpLt(aa, LFS(vul_ind), lo, hi);
        }
        
        Type max_vul = max(vul_age);
        for(int a=0;a<n_age;a++) vul(y,a,ff) = vul_age(a)/max_vul;
      } else { //if(vul_type(vul_ind) == -2) { // Free parameters
        for(int a=0;a<n_age;a++) vul(y,a,ff) = invlogit(vul_par(a, vul_ind));
      } //else { // Age-specific index - superseded by free parameters long ago
        //vul(y,vul_type(vul_ind)-1,ff) = 1;
      //}
    }
  }
  return vul;
}

template<class Type>
array<Type> calc_ivul(matrix<Type> vul_par, vector<int> vul_type, vector<Type> lbinmid, int n_y, int n_age,
                      vector<matrix<Type> > PLA, vector<Type> &LFS, vector<Type> &L5,
                      vector<Type> &Vmaxlen, Type Linf, matrix<Type> mat, array<Type> fleet_var, 
                      matrix<Type> &vul_len, Type &prior, matrix<int> est_vul) {
  
  int nlbin = lbinmid.size();
  array<Type> vul(n_y, n_age, vul_type.size()); // Corresponding age based selectivity
  vul.setZero();
  
  vector<Type> sls(vul_type.size());
  vector<Type> srs(vul_type.size());
  Type maxage = Type(n_age) - 1;

  for(int ff=0;ff<vul_type.size();ff++) {

    if(vul_type(ff) == -4) { // B
      for(int y=0;y<n_y;y++) {
        for(int a=0;a<n_age;a++) vul(y,a,ff) = 1;
      }
    } else if(vul_type(ff) == -3) { // SSB
      for(int y=0;y<n_y;y++) {
        for(int a=0;a<n_age;a++) vul(y,a,ff) = mat(y,a);
      }
    } else if(vul_type(ff) == -2) { // free parameters
      for(int a=0;a<n_age;a++) {
        Type v = invlogit(vul_par(a,ff));
        if(est_vul(a,ff)) prior -= dbeta_(v, Type(1.01), Type(1.01), true) + log(v - v * v);
        for(int y=0;y<n_y;y++) vul(y,a,ff) = v;
      }
    } else if(vul_type(ff) > 0) { // Index mirrored to fleet
      vul.col(ff) = fleet_var.col(vul_type(ff) - 1);
    } else { // Logistic or dome vul_type %in% c(0, -1, -5, -6)
      if(est_vul(1,ff)) prior -= dnorm_(vul_par(1,ff), Type(0), Type(3), true);
      
      Type ilogit_x = invlogit(vul_par(0,ff));
      if(est_vul(0,ff)) {
        Type jac = ilogit_x - ilogit_x * ilogit_x;
        prior -= dbeta_(ilogit_x, Type(1.01), Type(1.01), true) + log(jac);
        //prior -= dnorm_(vul_par(0,ff), Type(0), Type(3), true);
      }
      if (vul_type(ff) == 0 || vul_type(ff) == -1) {
        LFS(ff) = ilogit_x * Linf;
      } else {
        LFS(ff) = ilogit_x * maxage;
      }
      L5(ff) = LFS(ff) - exp(vul_par(1,ff));
      sls(ff) = (LFS(ff) - L5(ff))/pow(-log2(0.05), 0.5);
      
      if(vul_type(ff) == -1 || vul_type(ff) == -6) { // Logistic
        Vmaxlen(ff) = 1;
      } else { // Dome
        Vmaxlen(ff) = invlogit(vul_par(2,ff));
        
        if(vul_type(ff) == 0) {
          srs(ff) = (Linf - LFS(ff))/pow(-log2(Vmaxlen(ff)), 0.5);
        } else {
          srs(ff) = (maxage - LFS(ff))/pow(-log2(Vmaxlen(ff)), 0.5);
        }
        
        if(est_vul(2,ff)) {
          Type jac = Vmaxlen(ff) - Vmaxlen(ff) * Vmaxlen(ff);
          prior -= dbeta_(Vmaxlen(ff), Type(1.01), Type(1.01), true) + log(jac);
        }
      }
      if (vul_type(ff) == 0 || vul_type(ff) == -1) { // Calculate length-based sel and convert to age
        for(int j=0;j<nlbin;j++) {
          Type lo = pow(2, -((lbinmid(j) - LFS(ff))/sls(ff) * (lbinmid(j) - LFS(ff))/sls(ff)));
          Type hi;
          if(vul_type(ff) == -1) {
            hi = 1;
          } else {
            hi = pow(2, -((lbinmid(j) - LFS(ff))/srs(ff) * (lbinmid(j) - LFS(ff))/srs(ff)));
          }
          vul_len(j,ff) = CppAD::CondExpLt(lbinmid(j), LFS(ff), lo, hi);
        }
        for(int y=0;y<n_y;y++) { // Ensure max(vul) = 1
          vector<Type> vul_age(n_age);
          vul_age.setZero();
          for(int a=0;a<n_age;a++) {
            for(int j=0;j<nlbin;j++) vul_age(a) += PLA(y)(a,j) * vul_len(j,ff);
          }
          Type max_vul = max(vul_age);
          for(int a=0;a<n_age;a++) vul(y,a,ff) = vul_age(a)/max_vul;
        }
      } else { // Calculate age-based sel
        vector<Type> vul_age(n_age);
        
        for(int a=0;a<n_age;a++) {
          Type aa = Type(a);
          Type lo = pow(2, -((aa - LFS(ff))/sls(ff) * (aa - LFS(ff))/sls(ff)));
          Type hi;
          if(vul_type(ff) == -6) {
            hi = 1;
          } else {
            hi = pow(2, -((aa - LFS(ff))/srs(ff) * (aa - LFS(ff))/srs(ff)));
          }
          vul_age(a) = CppAD::CondExpLt(aa, LFS(ff), lo, hi);
        }
        Type max_vul = max(vul_age);
        
        for(int y=0;y<n_y;y++) {
          for(int a=0;a<n_age;a++) vul(y,a,ff) = vul_age(a)/max_vul;
        }
      }
    }
    
  }
  return vul;
}


  
// Multinomial likelihood
template<class Type>
Type comp_multinom(array<Type> obs, array<Type> pred, matrix<Type> N, matrix<Type> N_samp, int y, int n_bin, int ff) {
  vector<Type> p_pred(n_bin);
  vector<Type> N_obs(n_bin);
  for(int bb=0;bb<n_bin;bb++) {
    p_pred(bb) = pred(y,bb,ff)/N(y,ff);
    N_obs(bb) = obs(y,bb,ff) * N_samp(y,ff);
  }
  Type log_like = dmultinom_(N_obs, p_pred, true);
  return log_like;
}

// Dirichlet multinomial linear version
template<class Type>
Type comp_dirmult1(array<Type> obs, array<Type> pred, matrix<Type> N, matrix<Type> N_samp, Type theta, int y, int n_bin, int ff) {
  vector<Type> alpha(n_bin);
  vector<Type> N_obs(n_bin);
  for(int bb=0;bb<n_bin;bb++) {
    alpha(bb) = theta * N_samp(y,ff) * pred(y,bb,ff)/N(y,ff);
    N_obs(bb) = obs(y,bb,ff) * N_samp(y,ff);
  }
  Type log_like = ddirmnom_(N_obs, alpha, true);
  return log_like;
}


// Dirichlet multinomial saturating version
template<class Type>
Type comp_dirmult2(array<Type> obs, array<Type> pred, matrix<Type> N, matrix<Type> N_samp, Type beta, int y, int n_bin, int ff) {
  vector<Type> alpha(n_bin);
  vector<Type> N_obs(n_bin);
  for(int bb=0;bb<n_bin;bb++) {
    alpha(bb) = beta * pred(y,bb,ff)/N(y,ff);
    N_obs(bb) = obs(y,bb,ff) * N_samp(y,ff);
  }
  Type log_like = ddirmnom_(N_obs, alpha, true);
  return log_like;
}

template<class Type>
Type comp_lognorm(array<Type> obs, array<Type> pred, matrix<Type> N, int y, int n_bin, int ff) {
  Type log_like = 0;
  for(int bb=0;bb<n_bin;bb++) {
    Type p_pred = pred(y,bb,ff)/N(y,ff);
    Type p_obs = obs(y,bb,ff);
    log_like += dnorm_(log(p_obs), log(p_pred), pow(0.02/p_obs, 0.5), true);
  }
  return log_like;
}

template<class Type>
Type comp_mvlogistic(array<Type> obs, array<Type> pred, matrix<Type> N, int n_y, int n_bin, int ff, Type p_min = 1e-8) {
  
  matrix<Type> tmp(n_y, n_bin);
  Type tau2 = 0;
  Type sum_count = 0; // (A - 1) * T
  
  vector<Type> A(n_y); // Number of age classes per year
  A.setZero();
  
  vector<Type> sum_term(n_y); // Annual mean residual deviation: mean(log(o/p))
  sum_term.setZero();
  
  for(int y=0;y<n_y;y++) {
    
    Type accum_obs = 0;
    Type accum_pred = 0;
    
    for(int bb=0;bb<n_bin;bb++) {
      Type p_pred = pred(y,bb,ff)/N(y,ff);
      tmp(y,bb) = CppAD::CondExpGt(obs(y,bb,ff), p_min, log(obs(y,bb,ff)) - log(p_pred), Type(0)); // Residual
      
      accum_obs += CppAD::CondExpLe(obs(y,bb,ff), p_min, obs(y,bb,ff), Type(0));
      accum_pred += CppAD::CondExpLe(obs(y,bb,ff), p_min, p_pred, Type(0));
      A(y) += CppAD::CondExpGt(obs(y,bb,ff), p_min, Type(1), Type(0)); // Zero if no data in year y
      
      sum_term(y) += tmp(y,bb);
    }
    
    A(y) += CppAD::CondExpGt(accum_obs, Type(0), Type(1), Type(0)); // One if no data in year y
    sum_term(y) += CppAD::CondExpGt(accum_obs, Type(0), log(accum_obs) - log(accum_pred), Type(0));
    sum_term(y) /= A(y);
    
    // Calculate eta iff A(y) > 1
    for(int bb=0;bb<n_bin;bb++) {
      tau2 += CppAD::CondExpGt(obs(y,bb,ff), p_min, (tmp(y,bb) - sum_term(y)) * (tmp(y,bb) - sum_term(y)), Type(0));
    }
    
    tau2 += CppAD::CondExpGt(A(y), Type(1), 
                             (log(accum_obs) - log(accum_pred) - sum_term(y)) * (log(accum_obs) - log(accum_pred) - sum_term(y)), 
                              Type(0));
    
    sum_count += CppAD::CondExpGt(A(y), Type(1), A(y) - 1, Type(0));
  }
  tau2 /= sum_count;
  
  Type log_like = -0.5 * sum_count * log(tau2) - 0.5 * sum_count;
  
  return log_like;
}


}
