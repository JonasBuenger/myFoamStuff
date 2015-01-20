#include"fvCFD.H"
#include "blockMatrixTools.H"

namespace myFvc{

template<class Type, class vectorType>
void div(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<GeometricField<typename innerProduct<vector, Type>::type, fvPatchField, volMesh> > tdivVf = fvc::div(vf);
    GeometricField<typename innerProduct<vector, Type>::type, fvPatchField, volMesh>& divVf = tdivVf();

    const int nComponents = pow(3,Type::rank-1);
    for(int i=0; i<nComponents; i++){

        blockMatrixTools::blockAdd(dirI+i, -f * divVf.internalField().component(i), B);

    }

}

}
