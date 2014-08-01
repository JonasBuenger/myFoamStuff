#ifndef MYIMPLICITGRADSCHEMES_H
#define MYIMPLICITGRADSCHEMES_H

#include "fvCFD.H"

namespace Foam{


namespace myImplicitGradSchemes{

tmp<fvScalarMatrix> coefficients(volScalarField&, int numComponent);

}


}
#endif // MYIMPLICITGRADSCHEMES_H
