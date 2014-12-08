#ifndef SP_H
#define SP_H

#include"fvCFD.H"
#include"subMatrix.H"
#include"subBlock.H"

namespace myFvm {

subMatrix<subBlock<1,1>, subBlock<1,1> > Sp(volScalarField&);

subMatrix<subBlock<3,3>, subBlock<3,1> > Sp(volVectorField&);

}

#endif
