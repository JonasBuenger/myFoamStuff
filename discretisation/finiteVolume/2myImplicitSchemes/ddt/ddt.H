#ifndef DDT_H
#define DDT_H

#include"fvCFD.H"
#include"subMatrix.H"
#include"subBlock.H"

namespace myFvm {

subMatrix<subBlock<1,1>, subBlock<1,1> > ddt(volScalarField&);

subMatrix<subBlock<3,3>, subBlock<3,1> > ddt(volVectorField&);

}

#endif
