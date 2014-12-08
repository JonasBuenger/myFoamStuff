#include"div.H"

#include"fvCFD.H"
#include"fvMesh.H"

namespace myFvm{

subMatrix<subBlock<1, 3>, subBlock<1, 1> > div(volVectorField& vf){

    subMatrix<subBlock<1, 3>, subBlock<1, 1> > matrix(vf.mesh());

    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme(
           surfaceInterpolationScheme<scalar>::New(
                   vf.mesh(),
                   vf.mesh().schemesDict().interpolationScheme(vf.name())
           )
    );
    const surfaceInterpolationScheme<scalar>& interpolationScheme = tinterpScheme();

    tmp<surfaceScalarField> tweights = interpolationScheme.weights(mag(vf));
    const surfaceScalarField& weights = tweights;


    // reset diag and source
    //diag_   = vectorField(mesh_.nCells(), pTraits<vector>::zero);
    //source_ = scalarField(mesh_.nCells(), pTraits<scalar>::zero);

    //if(mesh_.moving()){  // upper lower only need to be recomputed if mesh moves (?)
        // CONTRIBUTION INTERNAL FIELD

        for(int i=0; i<vf.mesh().owner().size(); i++){
            int o = vf.mesh().owner()[i];
            int n = vf.mesh().neighbour()[i];
            scalar w = weights.internalField()[i];
            vector s = vf.mesh().Sf()[i];

            matrix.diag(o).element(0,0) += (s*w).x();
            matrix.diag(o).element(0,1) += (s*w).y();
            matrix.diag(o).element(0,2) += (s*w).z();

            matrix.diag(n).element(0,0) -= (s*(1-w)).x();
            matrix.diag(n).element(0,1) -= (s*(1-w)).y();
            matrix.diag(n).element(0,2) -= (s*(1-w)).z();

            matrix.upper(i).element(0,0) = (s*(1-w)).x();
            matrix.upper(i).element(0,1) = (s*(1-w)).y();
            matrix.upper(i).element(0,2) = (s*(1-w)).z();

            matrix.lower(i).element(0,0) = -(s*w).x();
            matrix.lower(i).element(0,1) = -(s*w).y();
            matrix.lower(i).element(0,2) = -(s*w).z();

        }
    //}

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<vectorField> tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<vectorField> tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const vectorField& ic = tic();              // internal coefficient
        const vectorField& bc = tbc();              // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        const vectorField& sn = patch.Sf();         // patch normals

        forAll(patchVolField, faceI){               // loop all faces

            const label c = patch.faceCells()[faceI];
            const vector v = cmptMultiply(ic[faceI],sn[faceI]);

            matrix.diag(c).element(0,0) += v.x();
            matrix.diag(c).element(0,1) += v.y();
            matrix.diag(c).element(0,2) += v.z();

            matrix.source(c).element(0,0) -= bc[faceI] & sn[faceI];
       }

    }

    return matrix;

}


}
