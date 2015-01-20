#include "fvCFD.H"
#include "blockMatrixTools.H"

namespace myFvc{

template<class Type, class vectorType>
void gradDiv(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<GeometricField<typename innerProduct<vector,Type>::type, fvPatchField, volMesh> > tdivVf = fvc::div(vf);
    GeometricField<typename innerProduct<vector,Type>::type, fvPatchField, volMesh>& divVf = tdivVf();

    tmp<GeometricField<Type, fvPatchField, volMesh> > tgradDivVf = fvc::grad(divVf);
    GeometricField<Type, fvPatchField, volMesh>& gradDivVf = tgradDivVf();

    const int nComponents = pow(3,Type::rank);
    for(int i=0; i<nComponents; i++){

        blockMatrixTools::blockAdd(dirI+i, -1.0 * f * gradDivVf.internalField().component(i),B);

    }

}

template<class vectorType>
void gradDiv(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<GeometricField<vector, fvPatchField, volMesh> > tdivVf = fvc::div(vf);
    GeometricField<vector, fvPatchField, volMesh>& divVf = tdivVf();

    tmp<GeometricField<tensor, fvPatchField, volMesh> > tgradDivVf = fvc::grad(divVf);
    GeometricField<tensor, fvPatchField, volMesh>& gradDivVf = tgradDivVf();

    const int nComponents = 9;
    for(int i=0; i<nComponents; i++){

        blockMatrixTools::blockAdd(dirI+i, -1.0 * f * gradDivVf.internalField().component(i),B);

    }

}

template<class vectorType>
void gradDivSymm(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<volVectorField> tdivVf = fvc::div(vf);
    volVectorField& divVf = tdivVf();

    tmp<volTensorField> tgradDivVf = fvc::grad(divVf);
    volTensorField& gradDivVf = tgradDivVf();
    volSymmTensorField gradDivVfSymm = symm(gradDivVf);

    const int nComponents = 6;
    for(int i=0; i<nComponents; i++){

        blockMatrixTools::blockAdd(dirI+i, -1.0 * f * gradDivVfSymm.internalField().component(i),B);

    }
/*
    blockMatrixTools::blockAdd(dirI+0, -1.0 * f * gradDivVf.internalField().component(0),B);
    blockMatrixTools::blockAdd(dirI+1, -0.5 * f *(gradDivVf.internalField().component(1)+gradDivVf.internalField().component(3)),B);
    blockMatrixTools::blockAdd(dirI+2, -0.5 * f *(gradDivVf.internalField().component(2)+gradDivVf.internalField().component(6)),B);
    blockMatrixTools::blockAdd(dirI+3, -1.0 * f * gradDivVf.internalField().component(4),B);
    blockMatrixTools::blockAdd(dirI+4, -0.5 * f *(gradDivVf.internalField().component(5)+gradDivVf.internalField().component(7)),B);
    blockMatrixTools::blockAdd(dirI+5, -1.0 * f * gradDivVf.internalField().component(8),B);
*/

}

}
