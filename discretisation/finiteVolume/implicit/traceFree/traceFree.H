#ifndef TRACEFREEIMPLICIT_H
#define TRACEFREEIMPLICIT_H

#include"fvCFD.H"

namespace myFvm{

template<class vectorType>
void traceFree(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int eqnNumber, const int dirJ);

}

#include "./traceFree.C"

#endif
