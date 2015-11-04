//
//  elliptic.h
//  GLAMER
//
//  Created by bmetcalf on 8/5/14.
//
//

#ifndef __GLAMER__elliptic__
#define __GLAMER__elliptic__

#include <iostream>
#include "standard.h"
#include "slsimlib.h"
//#include "lens.h"

/** \brief Class to convert isotropic LensHalos into elliptical halos.
 *
 *  This class changes the mass distribution into an elliptical shape 
 *  using the approach of T. Schramm, 1990 A&A, 231,19. The mass is 
 *  conserved.  An isotropic LensHalo needs to be constructed first and 
 *  then it is past to Elliptic.
 *
 */
class Elliptic{
public:
  Elliptic(
           LensHalo *isohalo    /// isotropic halo to be ellipticized
           ,PosType fratio      /// axis ratio of ellipse
           ,PosType inclination /// position angle of ellipse
           )
  : inclination(inclination),isohalo(isohalo){
    a = isohalo->get_Rsize()*sqrt(fratio);
    b = isohalo->get_Rsize()/sqrt(fratio);
  };
  
  /// The deflection angle in solar mass/Mpc
  void alpha(PosType x[],PosType alpha[]);
private:

  PosType a;
  PosType b;
  PosType inclination;
  LensHalo *isohalo;

  struct DALPHAXDM{
    DALPHAXDM(PosType lambda,PosType a2,PosType b2,PosType X[],LensHalo* isohalo)
    :lambda(lambda),a2(a2),b2(b2),x(X),isohalo(isohalo){};
    
    PosType operator()(PosType m);
  private:
    PosType lambda;
    PosType a2;
    PosType b2;
    PosType *x;
    LensHalo* isohalo;
  };
  struct DALPHAYDM{
    DALPHAYDM(PosType lambda,PosType a2,PosType b2,PosType X[],LensHalo* isohalo)
    :lambda(lambda),a2(a2),b2(b2),x(X),isohalo(isohalo){};
    
    PosType operator()(PosType m);
  private:
    PosType lambda;
    PosType a2;
    PosType b2;
    PosType *x;
    LensHalo* isohalo;
  };
  
};


#endif /* defined(__GLAMER__elliptic__) */
