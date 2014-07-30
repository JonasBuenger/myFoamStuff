#include "interpolateThetaAtBoundary.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"
#include "fvcGrad.H"
#include "fvm.H"

scalarField Foam::OFinterpolate(const fvPatch& patch, const objectRegistry& db, scalarField boundaryTheta){

    const volScalarField& Theta = db.lookupObject<volScalarField>("Theta");
    const vectorField ThetaGradient = fvc::grad(Theta,"leastSquares");
    const vectorField& patchDeltas = patch.delta();     // cell-centre to face-centre vector

    for(int i=0; i<patch.size(); i++){
        boundaryTheta[i] = Theta[patch.faceCells()[i]]
                            + ( ThetaGradient[i] & patchDeltas[i] );
    }

    return boundaryTheta;

}
