#ifndef DIVEXPLICIT_H
#define DIVEXPLICIT_H

#include"fvCFD.H"

namespace myFvc{

template<class Type, class vectorType>
void div(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "./div.C"

#endif
