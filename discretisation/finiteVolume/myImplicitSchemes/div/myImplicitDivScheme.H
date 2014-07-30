#ifndef MYIMPLICITDIVSCHEME_H
#define MYIMPLICITDIVSCHEME_H

#include "myImplicitScheme.H"
#include "fvCFD.H"

namespace Foam{


class myImplicitDivScheme : public myImplicitScheme{

public:
    myImplicitDivScheme(const volVectorField);
    tmp<fvScalarMatrix> coefficients(int numComponent);
};


}
#endif // MYIMPLICITDIVSCHEME_H
