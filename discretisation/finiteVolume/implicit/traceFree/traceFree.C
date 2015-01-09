#include"fvCFD.H"
#include "../useful/interpolationWeights.H"

namespace myFvm{

template<class vectorType>
void traceFree(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int eqnNumber, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    forAll(diag, i){
        for(int j=0; j<vectorType::nComponents; j++){
            diag[i](eqnNumber, j) = 0;
        }
        diag[i](eqnNumber,dirJ+0) = 1;
        diag[i](eqnNumber,dirJ+3) = 1;
        diag[i](eqnNumber,dirJ+5) = 1;
    }

    forAll(upper, i){
        for(int j=0; j<vectorType::nComponents; j++){
            upper[i](eqnNumber, j) = 0;
        }
    }

    forAll(lower, i){
        for(int j=0; j<vectorType::nComponents; j++){
            lower[i](eqnNumber, j) = 0;
        }
    }

    forAll(B, i){
        B[i](eqnNumber) = 0;
    }

}

}
