#include"fvCFD.H"
#include "blockMatrixTools.H"

namespace myFvc{

template<class vectorType>
void divDivI(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    tmp<volVectorField> tdivVf = fvc::div(vf);
    volVectorField& divVf = tdivVf();

    tmp<volScalarField> tdivDivVf = fvc::div(divVf);
    volScalarField& divDivVf = tdivDivVf();

    blockMatrixTools::blockAdd(dirI+0, -1.0 * f * divDivVf.internalField().component(0),B);
    blockMatrixTools::blockAdd(dirI+3, -1.0 * f * divDivVf.internalField().component(0),B);
    blockMatrixTools::blockAdd(dirI+5, -1.0 * f * divDivVf.internalField().component(0),B);

}

}
