
namespace ns_VPA {

using namespace ns_SCA;

// Iterative solver for F in VPA model -
// if F_ya = Z_ya C_ya / (exp(Z_ya) - 1) N_y+1,a+1, then:
// g(x) = log(Z) + log(C) - log(exp(Z) - 1) - log(N) - x
// where x = log(F), and:
// g'(x) = exp(x)/Z - exp(x) exp(Z)/(exp(Z) - 1) - 1
template<class Type>
Type VPA_F(Type logF, Type M, Type CAA, Type N_next) {
  Type Z = exp(logF) + M;
  return log(Z) + log(CAA) - log(N_next) - log(exp(Z) - 1) - logF;
}

// Derivative of VPA_F
template<class Type>
Type deriv_VPA_F(Type logF, Type M) {
  Type Z = exp(logF) + M;
  Type ans = exp(logF)/Z;
  ans -= exp(logF) * exp(Z)/(exp(Z) - 1);
  return ans - 1;
}

// Newton solver for log(F)
template<class Type>
Type Newton_VPA_F(Type F_start, Type M, Type CAA, Type N_next, int nloop) {
  Type logF = log(F_start);
  for(int i=0;i<nloop;i++) {
    Type tmp = VPA_F(logF, M, CAA, N_next)/deriv_VPA_F(logF, M);
    logF -= tmp;
  }
  return exp(logF);
}

// Iterative solver for F for A-1 and A (maximum age as a plus-group) in VPA model -
// Solve for x = log(F_y,A-1) such that:
// h(x) = N_y,A exp(-Z_y,A) + N_y,A-1 exp(-Z_y,A-1) - N_y+1,A = 0
// where:
// N_y,A-1 = Z_y,A-1 / F_y,A-1 / (1 - exp(-Z_y,A-1)) * C_y,A-1
// N_y,A = Z_y,A / F_y,A / (1 - exp(-Z_y,A)) * C_y,A where F_y,A = phi * F_y,A-1
template<class Type>
Type VPA_F_plus(Type logF, Type phi, Type M1, Type M2, Type C1, Type C2, Type N_next) {
  Type F1 = exp(logF);
  Type Z1 = F1 + M1;
  Type F2 = phi * F1;
  Type Z2 = F2 + M2;
  Type N1 = Z1/F1/(1 - exp(-Z1)) * C1;
  Type N2 = Z2/F2/(1 - exp(-Z2)) * C2;
  return N1 * exp(-Z1) + N2 * exp(-Z2) - N_next;
}

// Derivative of h(x)
template<class Type>
Type deriv_VPA_F_plus(Type logF, Type phi, Type M1, Type M2, Type C1, Type C2) {
  Type F1 = exp(logF);
  Type Z1 = F1 + M1;
  Type F2 = phi * F1;
  Type Z2 = F2 + M2;

  Type a = Z1/(1 - exp(-Z1));
  Type a_deriv = 1 - exp(-Z1) - Z1 * exp(-Z1);
  a_deriv *= F1;
  a_deriv /= (1 - exp(-Z1)) * (1 - exp(-Z1));

  Type b = C1/F1;
  Type b_deriv = -C1/F1;

  Type cc = Z2/(1 - exp(-Z2));
  Type cc_deriv = 1 - exp(-Z2) - Z2 * exp(-Z2);
  cc_deriv *= F2;
  cc_deriv /= (1 - exp(-Z2)) * (1 - exp(-Z2));

  Type d = C2/F2;
  Type d_deriv = -C2/F2;

  Type ans = a_deriv * b + a * b_deriv;
  ans += cc_deriv * d + cc * d_deriv;
  return ans;
}

// Newton solver for plus-group F
template<class Type>
Type Newton_VPA_F_plus(Type F_start, Type phi, Type M1, Type M2, Type CAA1, Type CAA2, Type N_next, int nloop) {
  Type logF = log(F_start);
  for(int i=0;i<nloop;i++) {
    Type tmp = VPA_F_plus(logF, phi, M1, M2, CAA1, CAA2, N_next);
    tmp /= deriv_VPA_F_plus(logF, phi, M1, M2, CAA1, CAA2);
    logF -= tmp;
  }
  return exp(logF);
}

}
