#ifndef LAPLACIAN_H
#define LAPLACIAN_H

#include"fvCFD.H"

namespace myFvm {

template<class vectorType>
void laplacian(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void laplacian(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void laplacian(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void laplacian(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void laplacian(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void laplacian(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

}

#include "laplacian.C"

#endif
