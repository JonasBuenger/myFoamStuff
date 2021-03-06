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

tmp<surfaceScalarField> Foam::myImplicitDivScheme::getWeights(){

    // interpolation scheme
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme(
           surfaceInterpolationScheme<scalar>::New(
                   mesh_,
                   mesh_.schemesDict().interpolationScheme(vf_.name())
           )
    );
    const surfaceInterpolationScheme<scalar>& interpolationScheme = tinterpScheme();

    tmp<surfaceScalarField> tweights = interpolationScheme.weights(mag(vf_));
    return tweights;

}

void
Foam::myImplicitDivScheme::updateCoeffs(){

    // interpolation weights
    tmp<surfaceScalarField> tweights = getWeights();
    const surfaceScalarField& weights = tweights;

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
            upper_[i] =  s*(1-w);
            lower_[i] = -s*w;

        }
    //}

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vf_boundFace = ic * vf_boundCell  +  bc

    vf_.boundaryField().updateCoeffs();
    forAll(vf_.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf_.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<vectorField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<vectorField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const vectorField& ic = tic();     // internal coefficient
        const vectorField& bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        const vectorField& sn = patch.Sf();   // patch normals

        forAll(patchVolField, faceI){ // loop all faces

            const label c = patch.faceCells()[faceI];

            diag_[c]   += cmptMultiply(ic[faceI],sn[faceI]);
            source_[c] -= bc[faceI] & sn[faceI];

        }

    }
    /*
    diag_ = vector(1,1,1);
    upper_ = vector(2,2,2);
    lower_ = vector(3,3,3);
    source_ = 4;
    */
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
