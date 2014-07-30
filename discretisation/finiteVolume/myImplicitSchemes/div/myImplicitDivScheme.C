#include "myImplicitDivScheme.H"

#include "fvCFD.H"
#include "fvMatrix.H"

Foam::myImplicitDivScheme::myImplicitDivScheme(const volVectorField vf)
    :myImplicitScheme(vf)
{
    Info << "bla" << endl;
};


tmp<fvScalarMatrix>
Foam::myImplicitDivScheme::coefficients(int numComponent){

    // FVMATRIX
    volScalarField vf0 = vf_.component(numComponent);

    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix
        (
            vf0,
            vf_.dimensions()/dimLength
        )
    );
    fvScalarMatrix& fvm = tfvm();

    // INTERPOLATION SCHEME
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
            surfaceInterpolationScheme<scalar>::New(
                    vf_.mesh(),
                    mesh_.schemesDict().interpolationScheme(vf_.name())
            )
    );

    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf_.component(numComponent));
    const surfaceScalarField& weights = tweights();

    // COEFFS AND SOURCE
    scalarField& diag = fvm.diag();
    scalarField& upper = fvm.upper();
    scalarField& lower = fvm.lower();
    scalarField& source = fvm.source();

    // CONTRIBUTION INTERNAL FIELD
    for(int i=0; i<mesh_.owner().size(); i++){
        int o = mesh_.owner()[i];
        int n = mesh_.neighbour()[i];
        scalar w = weights.internalField()[i];
        scalar s = mesh_.Sf()[i].component(numComponent);

        diag[o] +=  s*w;
        diag[n] -=  s*(1-w);
        lower[i] = -s*w;
        upper[i] =  s*(1-w);

    }

    // CONTRIBUTION BOUNDARY FIELD ...      (update coeffs in solver!!!)
    //      vf_patchFace = ic * vf_boundaryCell  +  bc
    // vf_.boundaryField().updateCoeffs();

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

            diag[c]   += ic[faceI].component(numComponent) * sn[faceI].component(numComponent);
            source[c] -= bc[faceI].component(numComponent) * sn[faceI].component(numComponent);

        }

    }

    return tfvm;
};
