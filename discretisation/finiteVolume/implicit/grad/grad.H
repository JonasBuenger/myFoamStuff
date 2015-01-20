#ifndef GRADIMPLICIT_H
#define GRADIMPLICIT_H

#include"fvCFD.H"

namespace myFvm{

template<class vectorType>
void gaussGrad(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

template<class vectorType>
void gaussSymmGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

template<class vectorType>
void gaussGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "./grad.C"

#endif
