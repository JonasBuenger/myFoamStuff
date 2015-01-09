
#include "blockVectorNMatrices.H"
#include"fvCFD.H"
#include "fvMesh.H"

namespace myFvm{

template<class vectorType>
void ddt(volScalarField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, int dir, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    const scalar rDeltaT = 1.0/vf.mesh().time().deltaT().value();
    const scalarField& vf0 = vf.oldTime().internalField();

    if (!updateOnlyRHS){

        forAll(diag, i){

                diag[i](dir, dir) += f * rDeltaT;

        }

    }

    forAll(B, i){

        B[i](dir) += f * rDeltaT * vf0[i];

    }

}

template<class Type, class vectorType>
void ddt(GeometricField<Type, fvPatchField, volMesh>& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
         const int dir, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    const scalar rDeltaT = 1.0/vf.mesh().time().deltaT().value();
    const Field<Type>& vf0 = vf.oldTime().internalField();

    int nComp = Type::nComponents;

    if (!updateOnlyRHS){

        forAll(diag, i){

            for(int j = 0; j<nComp; j++){

                diag[i](dir + j, dir + j) += f * rDeltaT;

            }

        }

    }

    forAll(B, i){

        for(int j = 0; j<nComp; j++){

            B[i](dir + j) += f * rDeltaT * vf0[i].component(j);

        }

    }

}


}
