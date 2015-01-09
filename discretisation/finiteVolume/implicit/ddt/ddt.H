#ifndef DDT_H
#define DDT_H

#include "fvCFD.H"
#include "blockVectorNMatrices.H"

namespace myFvm {

template<class Type, class vectorType>
void ddt(GeometricField<Type, fvPatchField, volMesh>& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
         const int dir, const bool updateOnlyRHS = false);

template<class vectorType>
void ddt(volScalarField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, int dir, const bool updateOnlyRHS = false);

}

#include "ddt.C"

#endif
