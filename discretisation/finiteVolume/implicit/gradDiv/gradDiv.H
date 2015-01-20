#ifndef GRADDIVIMPLICIT_H
#define GRADDIVIMPLICIT_H

#include"fvCFD.H"

namespace myFvm{

template<class vectorType>
void gaussGradDivSymm(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "./gradDiv.C"

#endif
