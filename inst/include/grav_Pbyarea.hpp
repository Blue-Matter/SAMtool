
#ifndef grav_Pbyarea_hpp
#define grav_Pbyarea_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type grav_Pbyarea(objective_function<Type> *obj) {
  using namespace density;

  DATA_VECTOR( fracs );
  DATA_VECTOR( prob );
  DATA_INTEGER( nareas );

  PARAMETER_VECTOR( log_visc );
  PARAMETER_VECTOR( log_grav );

  Type sigma = Type(0);
  sigma=0.1;

  // -- Declarations
  matrix<Type> grav(nareas,nareas);
  matrix<Type> mov(nareas,nareas);
  matrix<Type> transN(nareas,nareas);
  vector<Type> idist(nareas);
  vector<Type> tdist(nareas);
  vector<Type> gsum(nareas);
  vector<Type> Nsum(nareas);

  // Zero inits
  grav.setZero();
  mov.setZero();
  idist.setZero();
  tdist.setZero();
  gsum.setZero();
  Nsum.setZero();
  transN.setZero();

  // Map out gravity terms accounting for viscosity, gravity of area 1 is fixed to zero
  for(int af=0; af<nareas; af++){ // area from
    for(int at=1; at<nareas; at++){ // area to
      grav(af,at) = log_grav(at-1);
    }
  }
  for(int af=0; af<nareas; af++){
    grav(af,af)+=log_visc(af);
  }

  // Calculate logit fractions (movement to area from area)
  for(int af=0; af<nareas; af++){

    for(int at=0; at<nareas; at++){
      gsum(af)+=exp(grav(af,at)); // sum up gravity terms by from area
    }

    for(int at=0; at<nareas; at++){
      mov(af,at)=exp(grav(af,at))/gsum(af); // calculate logit mov probs by row (area from)
    }

  }

  //std::cout<<mov<<std::endl;
  // Run a convergence to a stable distribution
  for(int af=0; af<nareas; af++){
    idist(af)=1.0/nareas;
  }

  for(int tt =0; tt<50; tt++){
     //tdist=idist*mov;
     //idist=tdist;
    for(int af=0; af<nareas; af++){
      for(int at=0; at<nareas; at++){
        transN(af,at)=idist(af)*mov(af,at);
      }
    }

    Nsum.setZero();

    for(int at=0; at<nareas; at++){
      for(int af=0; af<nareas; af++){
        Nsum(at)+=transN(af,at);
      }
      idist(at)=Nsum(at);
    }

  }

  Type nll = Type(0);
  for(int aa=0; aa<nareas; aa++){
   // nll-=dnorm(log(idist(aa)), log(fracs(aa)), sigma, true);
   // nll-=dnorm(log(mov(aa,aa)),log(prob(aa)),sigma,true);

    nll-=dnorm(idist(aa), fracs(aa), sigma, true);
    nll-=dnorm(mov(aa,aa),prob(aa),sigma,true);
  }

  //-------REPORTING-------//
  ADREPORT( log_grav );
  ADREPORT( log_visc );
  REPORT( idist );
  REPORT( transN );
  REPORT( grav );
  REPORT( mov );
  REPORT( grav );
  REPORT( log_grav );
  REPORT( log_visc );
  REPORT( prob );
  REPORT( fracs );

  return nll;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
