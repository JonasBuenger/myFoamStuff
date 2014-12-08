#ifndef GRAD_H
#define GRAD_H

#include"fvCFD.H"
#include"subMatrix.H"
#include"subBlock.H"

namespace myFvm{

subMatrix<subBlock<3,1>, subBlock<3,1> > grad(volScalarField& vf);

}

#endif
