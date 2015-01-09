#ifndef GRAD_H
#define GRAD_H

#include"fvCFD.H"

namespace myFvm{

template<class vectorType>
void gaussGrad(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussSymmGrad(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussGrad(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussGrad(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussSymmGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

}

#include "grad.C"

#endif
