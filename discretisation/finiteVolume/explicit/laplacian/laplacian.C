#include "fvCFD.H"
#include "blockMatrixTools.H"

namespace myFvc{

template<class Type, class vectorType>
void laplacian(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<GeometricField<Type, fvPatchField, volMesh> > tlaplacianVf = fvc::laplacian(vf);
    GeometricField<Type, fvPatchField, volMesh>& laplacianVf = tlaplacianVf();

    int nComponents = pow(3,Type::rank);
    if(Type::typeName == "symmTensor"){
        nComponents = 6;
    }
    for(int i=0; i<nComponents; i++){

        blockMatrixTools::blockAdd(dirI+i, -f * laplacianVf.internalField().component(i), B);

    }

}

}
