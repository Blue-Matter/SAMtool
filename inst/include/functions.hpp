
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
vector<Type> calc_q(matrix<Type> I_y, vector<Type> B_y, matrix<Type> &Ipred, int nsurvey) {
  vector<Type> q(nsurvey);

  for(int sur=0;sur<nsurvey;sur++) {
    vector<Type> I_vec = I_y.col(sur);
    q(sur) = calc_q(I_vec, B_y);
    for(int y=0;y<I_y.rows();y++) Ipred(y,sur) = q(sur) * B_y(y);
  }
  return q;
}

template<class Type>
vector<Type> calc_q(matrix<Type> I_y, vector<Type> B_y, vector<Type> N_y, matrix<Type> &Ipred, int nsurvey,
                    vector<int> I_units) {
  vector<Type> q(nsurvey);
  for(int sur=0;sur<nsurvey;sur++) {
    vector<Type> I_vec = I_y.col(sur);
    if(I_units(sur)) {
      q(sur) = calc_q(I_vec, B_y);
    } else {
      q(sur) = calc_q(I_vec, N_y);
    }
    for(int y=0;y<I_y.rows();y++) {
      if(I_units(sur)) {
        Ipred(y,sur) = q(sur) * B_y(y);
      } else {
        Ipred(y,sur) = q(sur) * N_y(y);
      }
    }
  }
  return q;
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

#include "ns/ns_cDD.hpp"
#include "ns/ns_SCA.hpp"
#include "ns/ns_VPA.hpp"
#include "ns/ns_SRA_scope.hpp"
#include "ns/ns_SP.hpp"


// This is the Newton solver for the fleet-specific F in year y given the observed catches.
// Let vector x = log(F_y,f)
// Then g(x) = Cpred - Cobs = sum_a v_y,a,f F_y,f / Z_y,a * (1 - exp(-Z_y,a)) * N_y,a * wt_y,a - Cobs
// g'(x) = sum_a v * N * w * deriv where deriv is defined in the code below
// We iteratively solve for x where x_next = x_previous - g(x)/g'(x)
template<class Type>
vector<Type> Newton_SRA_F(matrix<Type> C_hist, matrix<Type> N, matrix<Type> M, matrix<Type> wt, matrix<Type> VB_out, array<Type> vul,
                          Type max_F, int y, int max_age, int nfleet, int nit_F, Type &penalty) {

  vector<Type> F_out(nfleet);
  vector<Type> x_loop(nfleet);
  vector<Type> C_hist_y = C_hist.row(y);
  for(int ff=0;ff<nfleet;ff++) { // Starting values of x = log(F)
    Type F_start = CppAD::CondExpGt(C_hist(y,ff)/VB_out(y,ff), Type(0.95), Type(3), -log(1 - C_hist(y,ff)/VB_out(y,ff)));
    x_loop(ff) = log(F_start);
  }

  for(int i=0;i<nit_F;i++) { // Loop for Newton-Raphson
    vector<Type> Z = M.row(y);
    matrix<Type> VB(max_age, nfleet);
    vector<Type> Cpred(nfleet);
    vector<Type> F_loop(nfleet);
    Cpred.setZero();

    for(int ff=0;ff<nfleet;ff++) {
      F_loop(ff) = exp(x_loop(ff));
      VB.col(ff) = N.row(y);
      for(int a=0;a<max_age;a++) {
        VB(a,ff) *= vul(y,a,ff) * wt(y,a);
        Z(a) += vul(y,a,ff) * F_loop(ff);
      }
    }

    for(int ff=0;ff<nfleet;ff++) {
      for(int a=0;a<max_age;a++) Cpred(ff) += VB(a,ff) * F_loop(ff) * (1 - exp(-Z(a)))/Z(a);
    }

    if(i<nit_F-1) {
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
#include "grav.hpp"
#include "grav_Pbyarea.hpp"
#include "SCA.hpp"
#include "SCA_Pope.hpp"
#include "SCA2.hpp"
#include "SP.hpp"
#include "SRA_scope.hpp"
#include "VPA.hpp"

