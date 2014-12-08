
#include "blockVectorNMatrices.H"
#include"fvCFD.H"
#include "fvMesh.H"

namespace myFvm{

template<class vectorType>
void ddt(volScalarField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, int dir){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    const scalar rDeltaT = 1.0/vf.mesh().time().deltaT().value();
    const scalarField& vf0 = vf.oldTime().internalField();

    forAll(diag, i){
        diag[i](dir, dir) += f * rDeltaT;
    }

    forAll(B, i){
        B[i](dir) += f * rDeltaT * vf0[i];
    }

}


template<class vectorType>
void ddt(volVectorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    const scalar rDeltaT = 1.0/vf.mesh().time().deltaT().value();
    const vectorField& vf0 = vf.oldTime().internalField();

    forAll(diag, i){
        diag[i](dir + 0, dir + 0) += f * rDeltaT;
        diag[i](dir + 1, dir + 1) += f * rDeltaT;
        diag[i](dir + 2, dir + 2) += f * rDeltaT;
    }

    forAll(B, i){
        B[i](dir + 0) += f * rDeltaT * vf0[i].component(0);
        B[i](dir + 1) += f * rDeltaT * vf0[i].component(1);
        B[i](dir + 2) += f * rDeltaT * vf0[i].component(2);
    }

}

template<class vectorType>
void ddt(volSymmTensorField& vf, const double f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();

    const scalar rDeltaT = 1.0/vf.mesh().time().deltaT().value();
    const symmTensorField& vf0 = vf.oldTime().internalField();

    forAll(diag, i){
        diag[i](dir + 0, dir + 0) += f * rDeltaT;
        diag[i](dir + 1, dir + 1) += f * rDeltaT;
        diag[i](dir + 2, dir + 2) += f * rDeltaT;
        diag[i](dir + 3, dir + 3) += f * rDeltaT;
        diag[i](dir + 4, dir + 4) += f * rDeltaT;
        diag[i](dir + 5, dir + 5) += f * rDeltaT;
    }

    forAll(B, i){
        B[i](dir + 0) += f * rDeltaT * vf0[i].component(symmTensor::XX);
        B[i](dir + 1) += f * rDeltaT * vf0[i].component(symmTensor::XY);
        B[i](dir + 2) += f * rDeltaT * vf0[i].component(symmTensor::XZ);
        B[i](dir + 3) += f * rDeltaT * vf0[i].component(symmTensor::YY);
        B[i](dir + 4) += f * rDeltaT * vf0[i].component(symmTensor::YZ);
        B[i](dir + 5) += f * rDeltaT * vf0[i].component(symmTensor::ZZ);
    }

}


template<class vectorType>
void ddt(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, int dir){

    ddt(vf, 1.0, M, B, dir);

}


template<class vectorType>
void ddt(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir){

    ddt(vf, 1.0, M, B, dir);

}

template<class vectorType>
void ddt(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir){

    ddt(vf, 1.0, M, B, dir);

}

template<class vectorType>
void ddtTraceFree(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dir){

    ddtTraceFree(vf, 1.0, M, B, dir);

}

}
