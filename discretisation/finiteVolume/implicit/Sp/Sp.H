#ifndef SPIMPLICIT_H
#define SPIMPLICIT_H

#include"fvCFD.H"

namespace myFvm {

template<class Type, class vectorType>
void Sp(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false);

}

#include "./Sp.C"

#endif
