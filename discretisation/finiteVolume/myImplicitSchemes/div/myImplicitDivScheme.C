#include "myImplicitDivScheme.H"

Foam::myImplicitDivScheme::myImplicitDivScheme(volVectorField& vf):
    vf_(vf),
    mesh_(vf.mesh())
{
   // initialize fields
   diag_   = vectorField(mesh_.nCells(), pTraits<vector>::zero);
   upper_  = vectorField(mesh_.owner().size(), pTraits<vector>::zero);
   lower_  = vectorField(mesh_.owner().size(), pTraits<vector>::zero);
   source_ = scalarField(mesh_.nCells(), pTraits<scalar>::zero);

   // update coefficients
   updateCoeffs();
}

void
Foam::myImplicitDivScheme::updateCoeffs(){

    // interpolation scheme
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme(
           surfaceInterpolationScheme<scalar>::New(
                   mesh_,
                   mesh_.schemesDict().interpolationScheme(vf_.name())
           )
    );
    const surfaceInterpolationScheme<scalar>& interpolationScheme = tinterpScheme();

    tmp<surfaceScalarField> tweights = interpolationScheme.weights(mag(vf_));
    const surfaceScalarField& weights = tweights();

    // reset diag and source
    diag_   = vectorField(mesh_.nCells(), pTraits<vector>::zero);
    source_ = scalarField(mesh_.nCells(), pTraits<scalar>::zero);

    //if(mesh_.moving()){  // upper lower only need to be recomputed if mesh moves (?)
        // CONTRIBUTION INTERNAL FIELD
        for(int i=0; i<mesh_.owner().size(); i++){
            int o = mesh_.owner()[i];
            int n = mesh_.neighbour()[i];
            scalar w = weights.internalField()[i];
            vector s = mesh_.Sf()[i];

            diag_[o] +=  s*w;
            diag_[n] -=  s*(1-w);
            lower_[i] = -s*w;
            upper_[i] =  s*(1-w);
        }
    //}

    // CONTRIBUTION BOUNDARY FIELD ...      (update coeffs in solver!!!)
    //
    //      vf_boundFace = ic * vf_boundCell  +  bc

    vf_.boundaryField().updateCoeffs();
    forAll(vf_.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf_.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<Field<vector> > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<Field<vector> > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const Field<vector> ic = tic();     // internal coefficient
        const Field<vector> bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        forAll(patchVolField, faceI){ // loop all faces

            label c = patch.faceCells()[faceI];

            diag_[c]   += Vector<scalar>(ic[faceI].x() * sn[faceI].x(),
                                         ic[faceI].y() * sn[faceI].y(),
                                         ic[faceI].z() * sn[faceI].z());
            source_[c] -= bc[faceI] & sn[faceI];

        }

    }

}

const vectorField&
Foam::myImplicitDivScheme::coeffsDiag(){
    return diag_;
}

const vectorField&
Foam::myImplicitDivScheme::coeffsUpper(){
    return upper_;
}

const vectorField&
Foam::myImplicitDivScheme::coeffsLower(){
    return lower_;
}

const scalarField&
Foam::myImplicitDivScheme::coeffsSource(){
    return source_;
}
