
#include"Sp.H"
#include"fvCFD.H"
#include"fvMesh.H"

namespace myFvm{

/*
template<class vectorType>
void Sp(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
        const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    Sp(vf,1,M,B,dirI,dirJ, updateOnlyRHS);

}

template<class vectorType>
void Sp(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
        const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    Sp(vf,1,M,B,dirI,dirJ, updateOnlyRHS);

}

template<class vectorType>
void Sp(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
        const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    Sp(vf,1,M,B,dirI,dirJ,updateOnlyRHS);

}


template<class vectorType>
void Sp(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
        const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    if(!updateOnlyRHS){

        // Identität s
        forAll(diag, i){
            diag[i](dirI+0, dirJ+0) += f;
            diag[i](dirI+1, dirJ+1) += f;
            diag[i](dirI+2, dirJ+2) += f;
        }

    }

}

template<class vectorType>
void Sp(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
        const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    if(!updateOnlyRHS){

        // Identität s
        forAll(diag, i){
            diag[i](dirI+0, dirJ+0) += f;
            diag[i](dirI+1, dirJ+1) += f;
            diag[i](dirI+2, dirJ+2) += f;
            diag[i](dirI+3, dirJ+3) += f;
            diag[i](dirI+4, dirJ+4) += f;
            diag[i](dirI+5, dirJ+5) += f;
        }

    }

}
*/

template<class Type, class vectorType>
void Sp(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B,
        const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    const int nComp = Type::nComponents;

    if(!updateOnlyRHS){

        // Identität s
        forAll(diag, i){

            for(int j=0; j<nComp; j++){

                diag[i](dirI+j, dirJ+j) += f;

            }

        }

    }

}

}
