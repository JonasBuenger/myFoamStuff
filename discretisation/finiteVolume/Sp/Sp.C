
#include"Sp.H"
#include"fvCFD.H"
#include"fvMesh.H"

namespace myFvm{

template<class vectorType>
void Sp(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    Sp(vf,1,M,B,dirI,dirJ);

}

template<class vectorType>
void Sp(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    Sp(vf,1,M,B,dirI,dirJ);

}

template<class vectorType>
void Sp(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    Sp(vf,1,M,B,dirI,dirJ);

}

template<class vectorType>
void Sp(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    // Identität s
    forAll(diag, i){
        diag[i](dirI, dirJ) += f;
    }

}

template<class vectorType>
void Sp(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    // Identität s
    forAll(diag, i){
        diag[i](dirI+0, dirJ+0) += f;
        diag[i](dirI+1, dirJ+1) += f;
        diag[i](dirI+2, dirJ+2) += f;
    }

}

template<class vectorType>
void Sp(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

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
