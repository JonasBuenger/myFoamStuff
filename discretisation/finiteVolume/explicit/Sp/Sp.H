#ifndef SP_H
#define SP_H

#include"fvCFD.H"

namespace myFvm {

template<class vectorType>
void Sp(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void Sp(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void Sp(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void Sp(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void Sp(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void Sp(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

}

#include "Sp.C"

#endif
