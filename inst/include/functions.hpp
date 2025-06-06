

// plogis
template<class Type>
Type plogis(Type x, Type location = 0, Type scale = 1) {
  return 1/(1 + exp(-(x - location)/scale));
}

template<class Type>
Type qlogis(Type p, Type location = 0, Type scale = 1) {
  return scale * log(p/(1 - p)) + location;
}



template<class Type>
Type invlogit2(Type x, Type ymin = 0, Type ymax = 1, Type y0 = 0.5, Type scale = 1) {
  Type location = log((ymax - ymin)/(y0 - ymin) - 1);
  location *= scale;
  Type res = plogis(x, location, scale);
  res *= ymax - ymin;
  res += ymin;
  return res;
}


template<class Type>
Type logit2(Type v, Type ymin = 0, Type ymax = 1, Type y0 = 0.5, Type scale = 1) {
  Type location = log((ymax - ymin)/(y0 - ymin) - 1);
  location *= scale;
  Type p = (v - ymin)/(ymax - ymin);
  return qlogis(p, location, scale);
}

//posfun from ADMB
template<class Type>
Type posfun(Type x, Type eps, Type &penalty) {
  Type denom = 2;
  denom -= x/eps;
  Type ans = CppAD::CondExpGe(x, eps, x, eps/denom);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * (x - eps) * (x - eps));
  return ans;
}

//posfun2 from ADMB
template<class Type>
Type posfun2(Type x, Type eps, Type &penalty) {
  Type x_eps = x/eps;
  Type denom = 1.;
  denom += x_eps;
  denom += pow(x_eps, 2);
  denom += pow(x_eps, 3);
  Type pen2 = pow(denom, -1);
  pen2 *= eps;
  Type ans = CppAD::CondExpGe(x, eps, x, pen2);
  penalty += CppAD::CondExpGe(x, eps, Type(0), 0.01 * pow(x - eps, 3));
  return ans;
}

// Dirichlet multinomial density function
template <class Type>
Type ddirmnom(vector<Type> x, vector<Type> alpha, int give_log=0) {
  Type logres = lgamma(alpha.sum()) + lgamma(x.sum()+1) - lgamma(alpha.sum() + x.sum());
  for(int i=0;i<x.size();i++) {
    logres += lgamma(x(i)+alpha(i)) - lgamma(alpha(i)) - lgamma(x(i)+1);
  }
  if(give_log) return logres;
  else return exp(logres);
}

// Shortened Dirichlet multinomial density function
template <class Type>
Type ddirmnom_(vector<Type> x, vector<Type> alpha, int give_log=0) {
  Type logres = lgamma(alpha.sum()) - lgamma(alpha.sum() + x.sum());
  for(int i=0;i<x.size();i++) {
    logres += lgamma(x(i)+alpha(i)) - lgamma(alpha(i));
  }
  if(give_log) return logres;
  else return exp(logres);
}

// Shortened multinomial density function
template <class Type>
Type dmultinom_(vector<Type> x, vector<Type> p, int give_log=0) {
  Type logres = (x*log(p)).sum();
  if(give_log) return logres;
  else return exp(logres);
}

// Shortened beta density function
template<class Type>
Type dbeta_(Type x, Type shape1, Type shape2, int give_log) {
  Type logres = log(x) * (shape1-1) + log(1-x) * (shape2-1);
  if(give_log) return logres;
  else return exp(logres);
}
VECTORIZE4_ttti(dbeta_)

// Shortened normal density function
template<class Type>
Type dnorm_(Type x, Type mean, Type sd, int give_log=0) {
  Type resid = (x - mean) / sd;
  Type logans = - Type(.5) * resid * resid;
  if(give_log) return logans; else return exp(logans);
}
VECTORIZE4_ttti(dnorm_)


// Calculates analytical solution of a lognormal variable
template<class Type>
Type calc_sigma(vector<Type> I_y, vector<Type> Ipred_y) {
  Type sum_square = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.size();y++) {
    if(!R_IsNA(asDouble(I_y(y))) && I_y(y)>0) {
      sum_square += (log(I_y(y)/Ipred_y(y))) * (log(I_y(y)/Ipred_y(y)));
      n_y += 1.;
    }
  }
  Type sigma = pow(sum_square/n_y, 0.5);
  return sigma;
}

// Calculates analytical solution of a lognormal variable
//template<class Type>
//Type calc_sigma(matrix<Type> I_y, matrix<Type> Ipred_y, int nsurvey) {
//  vector<Type> sigma(nsurvey);
//  for(int sur=0;sur<nsurvey;sur++) {
//    vector<Type> obs = I_y.col(sur);
//    vector<Type> pred = Ipred_y.col(sur);
//    sigma(sur) = calc_sigma(obs, pred);
//  }
//  return sigma;
//}

// Calculates analytical solution of catchability when conditioned on catch and
// index is lognormally distributed.
template<class Type>
Type calc_q(vector<Type> I_y, vector<Type> B_y) {
  Type num = 0.;
  Type n_y = 0.;

  for(int y=0;y<I_y.size();y++) {
    if(!R_IsNA(asDouble(I_y(y))) && I_y(y)>0) {
      num += log(I_y(y)/B_y(y));
      n_y += 1.;
    }
  }
  Type q = exp(num/n_y);
  return q;
}

template<class Type>
Type calc_q(vector<Type> I_y, vector<Type> B_y, vector<Type> &Ipred_y) {
  Type q = calc_q(I_y, B_y);
  for(int y=0;y<Ipred_y.size();y++) Ipred_y(y) = q * B_y(y);
  return q;
}


// For SP, SCA
template<class Type>
vector<Type> calc_q(matrix<Type> I_y, vector<Type> B_y, matrix<Type> &Ipred, int nsurvey) {
  vector<Type> q(nsurvey);

  for(int sur=0;sur<nsurvey;sur++) {
    vector<Type> I_vec = I_y.col(sur);
    q(sur) = calc_q(I_vec, B_y);
    for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q(sur) * B_y(y);
  }
  return q;
}

// For DD, cDD
template<class Type>
vector<Type> calc_q(matrix<Type> I_y, vector<Type> B_y, vector<Type> N_y, matrix<Type> &Ipred, int nsurvey,
                    vector<int> I_units, int n_y) {
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    vector<Type> I_vec = I_y.col(sur);
    if(I_units(sur)) {
      q(sur) = calc_q(I_vec, B_y);
    } else {
      q(sur) = calc_q(I_vec, N_y);
    }
    for(int y=0;y<n_y;y++) {
      if(I_units(sur)) {
        Ipred(y,sur) = q(sur) * B_y(y);
      } else {
        Ipred(y,sur) = q(sur) * N_y(y);
      }
    }
  }
  return q;
}


// For RCM, VPA
template<class Type>
Type calc_q(matrix<Type> I_y, matrix<Type> B_y, int sur, int ff, matrix<Type> &Ipred, vector<int> abs_I, int n_y) {
  Type q;
  if(abs_I(sur)) {
    q = Type(1);
  } else {
    vector<Type> I_vec = I_y.col(sur);
    vector<Type> B_vec = B_y.col(ff);
    q = calc_q(I_vec, B_vec);
  }
  for(int y=0;y<n_y;y++) Ipred(y,sur) = q * B_y(y,ff);
  return q;
}

// For RCM, DD, cDD
template<class Type>
Type calc_prior(matrix<int> use_prior, matrix<Type> prior_dist, Type R0x, Type h, int SR_type, Type log_M, vector<Type> q, Type rescale) {
  Type prior = 0;
  if(use_prior(0) == 1) { // Prior for R0 - normal on log_R0, log Jacobian transform = zero
    prior += dnorm_(R0x - log(rescale), prior_dist(0,0), prior_dist(0,1), true);
  } else if(use_prior(0) == 2) { // uniform on log-R0, log Jacobian transform = zero
    prior -= log(log(prior_dist(0,1)) - log(prior_dist(0,0)));
  } else if(use_prior(0) == 3) { // uniform on R0, log Jacobian transform = log(r) + log(x)
    prior += -log(prior_dist(0,1) - prior_dist(0,0)) - log(rescale) + log(R0x);
  }
  if(use_prior(1)) { // Prior for h
    if(SR_type) { // Beverton-Holt - beta on y = (h - 0.2)/0.8 with log Jacobian transform of inverse logit fn
      Type y = (h - 0.2)/0.8;
      prior += dbeta_(y, prior_dist(1,0), prior_dist(1,1), true) + log(y - y * y); 
    } else { // Ricker - normal on h with log Jacobian transform
      prior += dnorm_(h, prior_dist(1,0), prior_dist(1,1), true) + log(h - 0.2);
    }
  }
  if(use_prior(2)) { // Prior for constant M - normal on log_M
    prior += dnorm_(log_M, prior_dist(2,0), prior_dist(2,1), true);
  }
  for(int i=3;i<use_prior.size();i++) { // Prior for q - normal on log(q)
    if(use_prior(i)) prior += dnorm_(log(q(i-3)), prior_dist(i,0), prior_dist(i,1), true);
  }
  return prior;
}

//////////// Functions for cDD.h, DD.h, SCA.h
template<class Type>
Type BH_SR(Type SSB, Type h, Type R0, Type SSB0) {
  Type Rpred = 4 * h * R0 * SSB;
  Type den = SSB0 * (1-h);
  den += (5*h-1) * SSB;
  Rpred /= den;
  return Rpred;
}

template<class Type>
Type Ricker_SR(Type SSB, Type h, Type R0, Type SSB0) {
  Type SSBPR0 = SSB0/R0;
  Type expon = 1;
  expon -= SSB/SSB0;
  expon *= 1.25;

  Type Rpred = pow(5*h, expon);
  Rpred *= SSB;
  Rpred /= SSBPR0;
  return Rpred;
}

template<class Type>
Type MesnilRochet_SR(Type x, Type gamma, Type Rmax, Type Shinge, int Sp = 1) {
  Type c1 = Shinge * Shinge + 0.25 * gamma * gamma;
  Type K = pow(c1, 0.5);
  Type beta = Rmax/(Shinge + K);
  Type Rpred;
  
  if (Sp) { // x is the number of spawners
    
    Type c2 = (x - Shinge) * (x - Shinge) + 0.25 * gamma * gamma;
    Type c3 = x + K - pow(c2, 0.5);
    Rpred = beta * c3;
    
  } else { // x is the spawners per recruit
    
    Type Se = 2*K/x/beta; // Equilibrium SSB
    Se -= 2 * (Shinge + K);
    
    Type Se_denom = 1/x/x/beta/beta;
    Se_denom -= 2/x/beta;
    
    Se /= Se_denom;
    
    Rpred = CppAD::CondExpGt(1/x, 2 * beta, Type(0), Se/x);
  }
  
  return Rpred;
}

#include "ns/ns_cDD.hpp"
#include "ns/ns_DD.hpp"
#include "ns/ns_SCA.hpp"
#include "ns/ns_VPA.hpp"
#include "ns/ns_RCM.hpp"
#include "ns/ns_SP.hpp"

// This is the Newton solver for the fleet-specific F in year y given the observed catches.
// Let vector x = log(F_y,f)
// Then g(x) = Cpred - Cobs = sum_a v_y,a,f F_y,f / Z_y,a * (1 - exp(-Z_y,a)) * N_y,a * wt_y,a,f - Cobs
// g'(x) = sum_a v * N * w * deriv where deriv is defined in the code below
// We iteratively solve for x where x_next = x_previous - g(x)/g'(x)
template<class Type>
vector<Type> Newton_F(matrix<Type> C_hist, matrix<Type> N, matrix<Type> M, array<Type> C_wt, matrix<Type> VB_out, array<Type> vul,
                      Type max_F, int y, int max_age, int nfleet, int n_itF, Type &penalty) {

  vector<Type> F_out(nfleet);
  vector<Type> x_loop(nfleet);
  vector<Type> C_hist_y = C_hist.row(y);
  for(int ff=0;ff<nfleet;ff++) { // Starting values of x = log(F)
    Type F_start = CppAD::CondExpGt(C_hist(y,ff)/VB_out(y,ff), Type(0.95), Type(3), -log(1 - C_hist(y,ff)/VB_out(y,ff)));
    x_loop(ff) = log(F_start);
  }

  for(int i=0;i<n_itF;i++) { // Loop for Newton-Raphson
    vector<Type> Z = M.row(y);
    matrix<Type> VB(max_age, nfleet);
    vector<Type> Cpred(nfleet);
    vector<Type> F_loop(nfleet);
    Cpred.setZero();

    for(int ff=0;ff<nfleet;ff++) {
      F_loop(ff) = exp(x_loop(ff));
      VB.col(ff) = N.row(y);
      for(int a=0;a<max_age;a++) {
        VB(a,ff) *= vul(y,a,ff) * C_wt(y,a,ff);
        Z(a) += vul(y,a,ff) * F_loop(ff);
      }
    }

    for(int ff=0;ff<nfleet;ff++) {
      for(int a=0;a<max_age;a++) Cpred(ff) += VB(a,ff) * F_loop(ff) * (1 - exp(-Z(a)))/Z(a);
    }

    if(i<n_itF-1) {
      vector<Type> Newton_fn = Cpred - C_hist_y;
      vector<Type> Newton_gr(nfleet);
      Newton_gr.setZero();
      for(int ff=0;ff<nfleet;ff++) {
        for(int a=0;a<max_age;a++) {
          Type deriv = F_loop(ff) * Z(a);
          deriv -= vul(y,a,ff) * F_loop(ff) * F_loop(ff);
          deriv *= 1 - exp(-Z(a));
          deriv += vul(y,a,ff) * F_loop(ff) * F_loop(ff) * Z(a) * exp(-Z(a));
          deriv /= Z(a) * Z(a);
          Newton_gr(ff) += VB(a,ff) * deriv;
        }
      }
      x_loop -= Newton_fn/Newton_gr; // Newton-Raphson

    } else {
      for(int ff=0;ff<nfleet;ff++) {
        Type tmp = max_F - F_loop(ff);
        F_out(ff) = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), F_loop(ff)); //posfun
      }
    }
  }
  return F_out;
}



#include "cDD.hpp"
#include "DD.hpp"
#include "SCA.hpp"
#include "SP.hpp"
#include "RCM.hpp"
#include "VPA.hpp"

