#ifndef DIVDIVEXPLICIT_H
#define DIVDIVEXPLICIT_H

#include"fvCFD.H"

namespace myFvc{

template<class vectorType>
void divDivI(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "divDivI.C"

#endif
