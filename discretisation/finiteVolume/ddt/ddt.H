#ifndef DDT_H
#define DDT_H

#include"fvCFD.H"
#include "blockVectorNMatrices.H"

namespace myFvm {

template<class vectorType>
void ddt(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, int dir);

template<class vectorType>
void ddt(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir);

template<class vectorType>
void ddt(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir);

template<class vectorType>
void ddt(volScalarField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, int dir);

template<class vectorType>
void ddt(volVectorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir);

template<class vectorType>
void ddt(volSymmTensorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir);

}

#include "ddt.C"

#endif
