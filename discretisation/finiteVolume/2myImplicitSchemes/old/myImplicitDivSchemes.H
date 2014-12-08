#ifndef MYIMPLICITDIVSCHEMES_H
#define MYIMPLICITDIVSCHEMES_H

#include "fvCFD.H"

namespace Foam{


namespace myImplicitDivSchemes{

tmp<fvScalarMatrix> coefficients(const volVectorField&, int numComponent);

}


}
#endif // MYIMPLICITDIVSCHEMES_H
