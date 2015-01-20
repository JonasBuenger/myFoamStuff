#ifndef GRADDIVEXPLICIT_H
#define GRADDIVEXPLICIT_H

#include"fvCFD.H"

namespace myFvc{

template<class Type, class vectorType>
void gradDiv(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

template<class vectorType>
void gradDivSymm(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "./gradDiv.C"

#endif
