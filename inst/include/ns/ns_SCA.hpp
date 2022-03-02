
namespace ns_SCA {

template<class Type>
vector<Type> calc_NPR(Type F, vector<Type> vul, vector<Type> M, int n_age, int Pope = 0) {
  vector<Type> NPR(n_age);
  NPR(0) = 1.;
  if(Pope) {
    for(int a=1;a<n_age;a++) NPR(a) = NPR(a-1) * exp(-M(a-1)) * (1 - vul(a-1) * F);
    NPR(n_age-1) /= 1 - exp(-M(n_age-1)) * (1 - vul(n_age-1) * F); // Plus-group
  } else {
    for(int a=1;a<n_age;a++) NPR(a) = NPR(a-1) * exp(-vul(a-1) * F - M(a-1));
    NPR(n_age-1) /= 1 - exp(-vul(n_age-1) * F - M(n_age-1)); // Plus-group
  }
  return NPR;
}

template<class Type>
vector<Type> calc_NPR(Type F, vector<Type> vul, Type M, int n_age, int Pope = 0) {
  vector<Type> Mv(n_age);
  Mv.fill(M);
  return calc_NPR(F, vul, Mv, n_age, Pope);
}

template<class Type>
vector<Type> calc_NPR(Type F, vector<Type> vul, matrix<Type> M, int n_age, int y, int Pope = 0) {
  vector<Type> Mv(n_age);
  Mv = M.row(y);
  return calc_NPR(F, vul, Mv, n_age, Pope);
}

template<class Type>
Type sum_EPR(vector<Type> NPR, vector<Type> weight, vector<Type> mat) {
  Type answer = 0.;
  for(int a=0;a<NPR.size();a++) answer += NPR(a) * weight(a) * mat(a);
  return answer;
}

template<class Type>
Type sum_BPR(vector<Type> NPR, vector<Type> weight) {
  Type answer = 0.;
  for(int a=0;a<NPR.size();a++) answer += NPR(a) * weight(a);
  return answer;
}

template<class Type>
Type sum_VBPR(vector<Type> NPR, vector<Type> weight, vector<Type> vul) {
  Type answer = 0.;
  for(int a=0;a<NPR.size();a++) answer += NPR(a) * weight(a) * vul(a);
  return answer;
}

template<class Type>
vector<Type> calc_logistic_vul(vector<Type> vul_par, int n_age, Type &prior) {
  vector<Type> vul(n_age);
  Type max_age = n_age - 1;
  Type a_95 = invlogit(vul_par(0)) * 0.75 * max_age;
  Type a_50 = a_95 - exp(vul_par(1));

  prior -= dnorm_(vul_par(0), Type(0), Type(3), true) + dnorm_(vul_par(1), Type(0), Type(3), true);

  for(int a=0;a<n_age;a++) {
    vul(a) = 1/(1 + exp(-log(Type(19.0)) * (Type(a) - a_50)/(a_95 - a_50)));
  }
  vul /= max(vul);

  return vul;
}

template<class Type>
Type dnorm_vul(Type x, Type mu, Type sd) {
  Type res = -0.5;
  Type resid = (x-mu)/sd;
  res *= resid * resid;
  return exp(res);
}

template<class Type>
vector<Type> calc_dome_vul(vector<Type> vul_par, int n_age, Type &prior) {
  vector<Type> vul(n_age);

  Type max_age = n_age - 1;
  Type a_full = invlogit(vul_par(0)) * 0.75 * max_age;
  Type a_50 = a_full - exp(vul_par(1));
  Type a_full2 = invlogit(vul_par(2));
  a_full2 *= max_age - a_full;
  a_full2 += a_full;
  Type vul_max = invlogit(vul_par(3));

  prior -= dnorm_(vul_par(0), Type(0), Type(3), true) + dnorm_(vul_par(1), Type(0), Type(3), true);
  prior -= dbeta_(vul_max, Type(1.01), Type(1.01), true) + log(vul_max - vul_max * vul_max);

  Type var_asc =(a_50 - a_full) * (a_50 - a_full);
  var_asc /= log(Type(4));

  Type var_des = (max_age - a_full2) * (max_age - a_full2);
  var_des /= -2 * log(vul_max);

  Type sd_asc = pow(var_asc, 0.5);
  Type sd_des = pow(var_des, 0.5);

  for(int a=0;a<n_age;a++) {
    Type aa = a;
    Type vul_asc = dnorm_vul(aa, a_full, sd_asc);
    Type vul_des = dnorm_vul(aa, a_full2, sd_des);

    vul(a) = CppAD::CondExpLe(aa, a_full, vul_asc, CppAD::CondExpLe(aa, a_full2, Type(1), vul_des));
  }
  vul /= max(vul);

  return vul;
}

template<class Type>
Type dlnorm_comp(vector<Type> obs, vector<Type> pred) {
  Type log_lik = 0.;
  for(int a=0;a<obs.size();a++) log_lik += dnorm_(log(obs(a)), log(pred(a)), 0.1/pow(pred(a), 0.5), true);
  return log_lik;
}

template<class Type>
Type calc_M_eq(Type F_eq, Type B0, Type R0, vector<Type> M_bounds, vector<Type> vul, vector<Type> weight, 
               int n_age, int Pope) {
  Type D = 0.4;
  Type Meq;
  for(int i=0;i<20;i++) {
    Meq = CppAD::CondExpLe(D, Type(1), M_bounds(0) + (M_bounds(1) - M_bounds(0)) * (1 - D), M_bounds(0));
    vector<Type> NPR = calc_NPR(F_eq, vul, Meq, n_age, Pope);
    Type B_eq = R0 * sum_BPR(NPR, weight);
    D = B_eq/B0;
  }
  return Meq;
}



}
