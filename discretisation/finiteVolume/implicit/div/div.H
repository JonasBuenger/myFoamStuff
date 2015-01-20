#ifndef DIVIMPLICIT_H
#define DIVIMPLICIT_H

#include"fvCFD.H"

namespace myFvm{

//template<class vectorType>
//void gaussDiv(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

//template<class vectorType>
//void gaussDiv(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

//template<class vectorType>
//void gaussDiv(volVectorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

//template<class vectorType>
//void gaussDiv(volSymmTensorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

//template<class vectorType>
//void gaussDiv(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

template<class Type, class vectorType>
void gaussDiv(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "./div.C"

#endif
