#ifndef DIV_H
#define DIV_H

#include"fvCFD.H"

namespace myFvm{

template<class vectorType>
void gaussDiv(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussDiv(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussDiv(volVectorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

template<class vectorType>
void gaussDiv(volSymmTensorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ);

}

#include"div.C"

#endif
