
namespace ns_SP {

// Sub-annual time steps for SP and iteratively find the F that predicts the observed catch
template<class Type>
Type SP_F(Type U_start, Type C_hist, Type MSY, Type K, Type n, Type nterm, Type dt, int nstep, int n_itF,
          vector<Type> &Cpred, vector<Type> &B, int y, Type &penalty) {
  Type F;
  Type F_out;
  Type B_out;
  if(nstep > 1) {
    F = -log(1 - U_start);
    for(int i=0;i<n_itF;i++) {
      Type Catch = 0;
      Type B_next = B(y);

      for(int seas=0;seas<nstep;seas++) {
        Type SP = CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B_next / K * log(B_next/K),
                                   nterm/(n-1) * MSY * (B_next/K - pow(B_next/K, n))) - F * B_next;
        SP *= dt;
        Catch += F * B_next * dt;
        B_next += SP;
      }
      F *= C_hist/Catch;
    }
    F_out = CppAD::CondExpLt(3 - F, Type(0), 3 - posfun(3 - F, Type(0), penalty), F);
    
    Type Catch = 0;
    Type B_next = B(y);
    for(int seas=0;seas<nstep;seas++) {
      Type SP = CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B_next / K * log(B_next/K),
                                 nterm/(n-1) * MSY * (B_next/K - pow(B_next/K, n))) - F_out * B_next;
      SP *= dt;
      Catch += F_out * B_next * dt;
      B_next += SP;
    }
    B_out = B_next;
    Cpred(y) = Catch;
  } else {
    F = C_hist/B(y);
    F_out = CppAD::CondExpLt(1 - F, Type(0.025), 1 - posfun(1 - F, Type(0.025), penalty), F);
    Cpred(y) = F_out * B(y);
    B_out = B(y) + CppAD::CondExpEq(n, Type(1), -exp(Type(1.0)) * MSY * B(y) / K * log(B(y)/K),
              nterm/(n-1) * MSY * (B(y)/K - pow(B(y)/K, n))) - F_out * B(y);
  }
  B(y+1) = CppAD::CondExpGt(B_out, Type(1e-8), B_out, posfun(B_out, Type(1e-8), penalty));
  return F_out;
}



}
