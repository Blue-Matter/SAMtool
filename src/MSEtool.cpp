#define TMB_LIB_INIT R_init_MSEtool
#include <TMB.hpp>
#include "../inst/include/functions.hpp"

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_STRING(model);

  if(model == "DD") {
    return DD(this);
  } else if(model =="SP") {
    return SP(this);
  } else if(model == "grav") {
    return grav(this);
  } else if(model == "grav_Pbyarea") {
	return grav_Pbyarea(this);
  } else if(model == "SCA") {
    return SCA(this);
  } else if(model == "SCA2") {
    return SCA2(this);
  } else if(model == "SCA_Pope") {
    return SCA_Pope(this);
  } else if(model == "VPA") {
    return VPA(this);
  } else if(model == "cDD") {
    return cDD(this);
  } else if(model == "SRA_scope") {
    return SRA_scope(this);
  } else {
    error("No model found.");
  }

  return 0;
}
