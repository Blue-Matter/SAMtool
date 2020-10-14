
namespace ns_cDD {

template<class Type>
Type cDD_BPR(Type F, Type M, Type wk, Type Kappa, Type Winf) {
  Type Z = F + M;
  Type BPR = Kappa * Winf/Z + wk;
  BPR /= Z + Kappa;
  return BPR;
}


template<class Type>
Type cDD_R(Type BPR, Type Arec, Type Brec, int SR_type) {
  Type Req;
  if(SR_type) {
    Req = Arec * BPR - 1;
  } else {
    Req = log(Arec * BPR);
  }
  return Req/Brec/BPR;
}

template<class Type>
Type cDD_F(Type F_start, Type C_hist, Type M, Type Winf, Type Kappa, Type wk, vector<Type> &N, vector<Type> &B, vector<Type> &Cpred,
           vector<Type> &BPRinf, vector<Type> &Binf, vector<Type> &R, vector<Type> &Ninf, int nitF, int tt) {
  Type F = F_start;

  for(int i=0;i<nitF;i++) {
    Type Z_tt = F + M;
    BPRinf(tt) = cDD_BPR(F, M, wk, Kappa, Winf);
    Binf(tt) = BPRinf(tt) * R(tt);
    Ninf(tt) = R(tt)/Z_tt;

    Type Catch = B(tt) - Binf(tt) - (N(tt) - Ninf(tt)) * Kappa * Winf / (Z_tt + Kappa);
    Catch *= 1 - exp(-Z_tt - Kappa);
    Catch /= Z_tt + Kappa;
    Catch += Binf(tt) + (N(tt) - Ninf(tt)) * Kappa * Winf / (Z_tt + Kappa);
    Catch *= F;

    if(i==nitF-1) {
      Cpred(tt) = Catch;
    } else {
      F *= C_hist/Catch;
    }
  }

  return F;
}


}
