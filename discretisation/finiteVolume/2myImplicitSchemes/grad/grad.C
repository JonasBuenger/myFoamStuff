#include"grad.H"

#include"fvCFD.H"
#include"fvMesh.H"

namespace myFvm{

subMatrix<subBlock<3, 1>, subBlock<3, 1> > grad(volScalarField& vf){

    subMatrix<subBlock<3, 1>, subBlock<3, 1> > matrix(vf.mesh());

    // interpolation scheme
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme(
           surfaceInterpolationScheme<scalar>::New(
                   vf.mesh(),
                   vf.mesh().schemesDict().interpolationScheme(vf.name())
           )
    );
    const surfaceInterpolationScheme<scalar>& interpolationScheme = tinterpScheme();

    tmp<surfaceScalarField> tweights = interpolationScheme.weights(vf);
    const surfaceScalarField& weights = tweights();

    // reset diag and source
    //diag_   = vectorField(vf.mesh().nCells(), pTraits<vector>::zero);
    //source_ = vectorField(vf.mesh().nCells(), pTraits<vector>::zero);

    // contribution of internal Field ...
    for(int i=0; i<vf.mesh().owner().size(); i++)
    {
        int o = vf.mesh().owner()[i];
        int n = vf.mesh().neighbour()[i];
        scalar w = weights.internalField()[i];
        vector s = vf.mesh().Sf()[i];

        matrix.diag(o).element(0,0) += (s*w).x();
        matrix.diag(o).element(1,0) += (s*w).y();
        matrix.diag(o).element(2,0) += (s*w).z();

        matrix.diag(n).element(0,0) -= (s*(1-w)).x();
        matrix.diag(n).element(1,0) -= (s*(1-w)).y();
        matrix.diag(n).element(2,0) -= (s*(1-w)).z();

        matrix.lower(i).element(0,0) = (-s*w).x();
        matrix.lower(i).element(1,0) = (-s*w).y();
        matrix.lower(i).element(2,0) = (-s*w).z();

        matrix.upper(i).element(0,0) = (s*(1-w)).x();
        matrix.upper(i).element(1,0) = (s*(1-w)).y();
        matrix.upper(i).element(2,0) = (s*(1-w)).z();

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<scalar>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

        tmp<scalarField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<scalarField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const scalarField& ic = tic();     // internal coefficient
        const scalarField& bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        forAll(patchVolField, faceI){   // loop all faces

            label c = patch.faceCells()[faceI];     // boundary cell

            matrix.diag(c).element(0,0) += (ic[faceI]*sn[faceI]).x();
            matrix.diag(c).element(1,0) += (ic[faceI]*sn[faceI]).y();
            matrix.diag(c).element(2,0) += (ic[faceI]*sn[faceI]).z();

            matrix.source(c).element(0,0) -= (bc[faceI]*sn[faceI]).x();
            matrix.source(c).element(1,0) -= (bc[faceI]*sn[faceI]).y();
            matrix.source(c).element(2,0) -= (bc[faceI]*sn[faceI]).z();

        }

    }

    return matrix;

}


}

