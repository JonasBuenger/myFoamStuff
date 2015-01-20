#include"fvCFD.H"
#include "blockMatrixTools.H"

namespace myFvc{

template<class Type, class vectorType>
void grad(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh> > tgradVf = fvc::grad(vf);
    GeometricField<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>& gradVf = tgradVf();

    int nComponents = pow(3,Type::rank + 1);
    for(int i=0; i<nComponents; i++){

        blockMatrixTools::blockAdd(dirI+i, -f * gradVf.internalField().component(i), B);

    }

}

}
