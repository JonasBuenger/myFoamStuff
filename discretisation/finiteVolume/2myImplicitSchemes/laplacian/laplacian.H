#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include"fvCFD.H"
#include"subMatrix.H"
#include"subBlock.H"

namespace myFvm {

subMatrix<subBlock<1,1>, subBlock<1,1> > laplacian(volScalarField&);

subMatrix<subBlock<3,3>, subBlock<3,1> > laplacian(volVectorField&);

}

#endif
