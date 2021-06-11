
namespace ns_DD {

template<class Type>
Type Newton_F(vector<Type> C_hist, Type M, vector<Type> B, Type max_F, int y, int n_itF, Type &penalty) {
  
  Type F_out;
  Type F_start = CppAD::CondExpGt(C_hist(y)/B(y), Type(0.95), Type(3), -log(1 - C_hist(y)/B(y)));
  Type x_loop = log(F_start);
  
  for(int i=0;i<n_itF;i++) { // Loop for Newton-Raphson
    Type F_loop = exp(x_loop);
    Type Z = F_loop + M;
    Type Cpred = B(y) * F_loop * (1 - exp(-Z))/Z;
    
    if(i<n_itF-1) {
      Type Newton_fn = Cpred - C_hist(y);
      Type Newton_gr = F_loop * Z;
      Newton_gr -= F_loop * F_loop;
      Newton_gr *= 1 - exp(-Z);
      Newton_gr += F_loop * F_loop * Z * exp(-Z);
      Newton_gr /= Z * Z;
      Newton_gr *= B(y);
      
      x_loop -= Newton_fn/Newton_gr;
    } else {
      Type tmp = max_F - F_loop;
      F_out = CppAD::CondExpLt(tmp, Type(0), max_F - posfun(tmp, Type(0), penalty), F_loop); //posfun
    }
  }
  return F_out;
}


}
