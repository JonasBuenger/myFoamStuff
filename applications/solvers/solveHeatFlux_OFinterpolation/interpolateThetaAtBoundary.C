#include "interpolateThetaAtBoundary.H"
#include "interpolationCellPoint.H"
#include "volPointInterpolation.H"

scalarField Foam::OFinterpolate(const fvPatch& patch, const objectRegistry& db, scalarField boundaryTheta, const vectorField& normals){


    // KOMMENTAR: Aus irgendeinem Grund stimmen normals in dieser Funktion zugegriffen nicht mit dem übergebenen
    //            normals überein. Die in der Funktion über patch.nf() stimmen nicht mit normals aus updateCoeffs überein.

    const fvMesh& mesh = patch.boundaryMesh().mesh();
    interpolationCellPoint<scalar> interpObject(db.lookupObject<volScalarField>("Theta"));
    const vectorField& patchCenters = patch.Cf();

    for(int i=0; i<patch.Cf().size(); i++){

        const scalar d = 1.0/patch.deltaCoeffs()[i];

        const vector pt0 = patchCenters[i];
        const vector pt1 = pt0 - normals[i]*d;
        const vector pt2 = pt1 - normals[i]*d;

        const label cellPt1 = mesh.findCell(pt1);
        const label cellPt2 = mesh.findCell(pt2);

        const scalar valPt1 = interpObject.interpolate(pt1, cellPt1);
        const scalar valPt2 = interpObject.interpolate(pt2, cellPt2);

        // boundaryTheta[i] = valPt1 - d * (valPt1 - valPt2) / d;
        boundaryTheta[i] = valPt1 - (valPt2 - valPt1);

    }

    return boundaryTheta;

}
