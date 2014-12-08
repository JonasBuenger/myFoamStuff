#ifndef DIV_H
#define DIV_H

#include"fvCFD.H"
#include"subMatrix.H"
#include"subBlock.H"

namespace myFvm{

subMatrix<subBlock<1,3>, subBlock<1,1> > div(volVectorField& vf);

}

//#include"div.C"

#endif
