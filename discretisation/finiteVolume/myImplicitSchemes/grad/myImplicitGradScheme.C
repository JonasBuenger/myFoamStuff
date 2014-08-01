#include "myImplicitGradScheme.H"

Foam::myImplicitGradScheme::myImplicitGradScheme(volScalarField& vf):
    vf_(vf),
    mesh_(vf.mesh())
{
   // initialize fields
   diag_   = vectorField(mesh_.nCells(), pTraits<vector>::zero);
   upper_  = vectorField(mesh_.owner().size(), pTraits<vector>::zero);
   lower_  = vectorField(mesh_.owner().size(), pTraits<vector>::zero);
   source_ = vectorField(mesh_.nCells(), pTraits<vector>::zero);

   // update coefficients
   updateCoeffs();
}

void
Foam::myImplicitGradScheme::updateCoeffs(){

    // interpolation scheme
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme(
           surfaceInterpolationScheme<scalar>::New(
                   mesh_,
                   mesh_.schemesDict().interpolationScheme("grad(" + vf_.name() + ")")
           )
    );
    const surfaceInterpolationScheme<scalar>& interpolationScheme = tinterpScheme();

    tmp<surfaceScalarField> tweights = interpolationScheme.weights(mag(vf_));
    const surfaceScalarField& weights = tweights();

    // reset diag and source
    diag_   = vectorField(mesh_.nCells(), pTraits<vector>::zero);
    source_ = vectorField(mesh_.nCells(), pTraits<vector>::zero);

    // contribution of internal Field ...
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

    // CONTRIBUTION BOUNDARY FIELD ...      (update coeffs in solver!!!)
    //
    //      vf_boundFace = ic * vf_boundCell  +  bc

    vf_.boundaryField().updateCoeffs();
    forAll(vf_.boundaryField(), patchI){ // loop all patches

        const fvPatchField<scalar>& patchVolField = vf_.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<Field<scalar> > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<Field<scalar> > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const Field<scalar> ic = tic();     // internal coefficient
        const Field<scalar> bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        forAll(patchVolField, faceI){ // loop all faces

            label c = patch.faceCells()[faceI];     // boundary cell

            diag_[c]   += ic[faceI]*sn[faceI];
            source_[c] -= bc[faceI]*sn[faceI];
        }

    }

}


const vectorField&
Foam::myImplicitGradScheme::coeffsDiag(){
    return diag_;
}

const vectorField&
Foam::myImplicitGradScheme::coeffsUpper(){
    return upper_;
}

const vectorField&
Foam::myImplicitGradScheme::coeffsLower(){
    return lower_;
}

const vectorField&
Foam::myImplicitGradScheme::coeffsSource(){
    return source_;
}

