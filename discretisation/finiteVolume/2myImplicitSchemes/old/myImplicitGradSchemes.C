#include "myImplicitGradSchemes.H"


tmp<fvScalarMatrix>
Foam::myImplicitGradSchemes::coefficients(volScalarField& vf, int numComponent){

    Info << "Foam::myImplicitDivSchemes::coefficients" << endl;

    const fvMesh& mesh = vf.mesh();

    // FV-MATRIX
    tmp<fvScalarMatrix> tfvm
    (
        new fvScalarMatrix
        (
            vf,
            vf.dimensions()
        )
    );
    fvScalarMatrix& fvm = tfvm();


    // INTERPOLATION SCHEME
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
            surfaceInterpolationScheme<scalar>::New(
                    vf.mesh(),
                    mesh.schemesDict().interpolationScheme(vf.name())
            )
    );

    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    // COEFFS AND SOURCE
    scalarField& diag = fvm.diag();
    scalarField& upper = fvm.upper();
    scalarField& lower = fvm.lower();
    scalarField& source = fvm.source();

    // CONTRIBUTION INTERNAL FIELD
    for(int i=0; i<mesh.owner().size(); i++){
        int o = mesh.owner()[i];
        int n = mesh.neighbour()[i];
        scalar w = weights.internalField()[i];
        scalar s = mesh.Sf()[i].component(numComponent);
        scalar invV = 1.0 / mesh.V()[i];

        diag[o] +=  s *   w   * invV;
        diag[n] -=  s * (1-w) * invV;
        lower[i] = -s *   w   * invV;
        upper[i] =  s * (1-w) * invV;

    }

    // CONTRIBUTION BOUNDARY FIELD ...      (update coeffs in solver!!!)
    //      vfpatchFace = ic * vfboundaryCell  +  bc
    // vf.boundaryField().updateCoeffs();

    forAll(vf.boundaryField(), patchI){    // loop all patches

        const fvPatchField<scalar>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<Field<scalar> > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<Field<scalar> > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const Field<scalar> ic = tic();     // internal coefficient
        const Field<scalar> bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        forAll(patchVolField, faceI){   // loop all faces

            label c = patch.faceCells()[faceI];

            diag[c]   += ic[faceI] * sn[faceI].component(numComponent);
            source[c] -= bc[faceI] * sn[faceI].component(numComponent);

        }

    }

    return tfvm;
}


tmp<scalarField>
Foam::myImplicitGradSchemes::upper(volScalarField& vf, int numComponent){

    Info << "Foam::myImplicitDivSchemes::upper" << endl;

    const fvMesh& mesh = vf.mesh();
    tmp<scalarField> tupper = tmp<scalarField>(new scalarField(vf.size(), pTraits<scalar>::zero));
    scalarField& upper = tupper();


    // INTERPOLATION SCHEME
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
            surfaceInterpolationScheme<scalar>::New(
                    vf.mesh(),
                    mesh.schemesDict().interpolationScheme(vf.name())
            )
    );

    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    // CONTRIBUTION INTERNAL FIELD
    for(int i=0; i<mesh.owner().size(); i++){
        scalar w = weights.internalField()[i];
        scalar s = mesh.Sf()[i].component(numComponent);
        scalar invV = 1.0 / mesh.V()[i];

        upper[i] =  s * (1-w) * invV;
    }

    return tupper;
}


tmp<scalarField>
Foam::myImplicitGradSchemes::lower(volScalarField& vf, int numComponent){

    Info << "Foam::myImplicitDivSchemes::lower" << endl;

    const fvMesh& mesh = vf.mesh();
    tmp<scalarField> tlower = tmp<scalarField>(new scalarField(vf.size(), pTraits<scalar>::zero));
    scalarField& lower = tlower();


    // INTERPOLATION SCHEME
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
            surfaceInterpolationScheme<scalar>::New(
                    vf.mesh(),
                    mesh.schemesDict().interpolationScheme(vf.name())
            )
    );

    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    // CONTRIBUTION INTERNAL FIELD
    for(int i=0; i<mesh.owner().size(); i++){
        scalar w = weights.internalField()[i];
        scalar s = mesh.Sf()[i].component(numComponent);
        scalar invV = 1.0 / mesh.V()[i];

        lower[i] = -s *   w   * invV;
    }

    return tlower;
}



tmp<scalarField>
Foam::myImplicitGradSchemes::diag(volScalarField& vf, int numComponent){

    Info << "Foam::myImplicitDivSchemes::diag" << endl;

    const fvMesh& mesh = vf.mesh();

    tmp<scalarField> tdiag = tmp<scalarField>(new scalarField(vf.size(), pTraits<scalar>::zero));
    scalarField& diag = tdiag();

    // INTERPOLATION SCHEME
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
            surfaceInterpolationScheme<scalar>::New(
                    vf.mesh(),
                    mesh.schemesDict().interpolationScheme(vf.name())
            )
    );

    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf);
    const surfaceScalarField& weights = tweights();

    // CONTRIBUTION INTERNAL FIELD
    for(int i=0; i<mesh.owner().size(); i++){
        int o = mesh.owner()[i];
        int n = mesh.neighbour()[i];
        scalar w = weights.internalField()[i];
        scalar s = mesh.Sf()[i].component(numComponent);
        scalar invV = 1.0 / mesh.V()[i];

        diag[o] +=  s *   w   * invV;
        diag[n] -=  s * (1-w) * invV;
    }

    // CONTRIBUTION BOUNDARY FIELD ...      (update coeffs in solver!!!)
    //      vfpatchFace = ic * vfboundaryCell  +  bc
    // vf.boundaryField().updateCoeffs();

    forAll(vf.boundaryField(), patchI){    // loop all patches

        const fvPatchField<scalar>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<Field<scalar> > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<Field<scalar> > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const Field<scalar> ic = tic();     // internal coefficient
        const Field<scalar> bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        forAll(patchVolField, faceI){   // loop all faces

            label c = patch.faceCells()[faceI];
            diag[c]   += ic[faceI] * sn[faceI].component(numComponent);

        }

    }

    return tdiag;
}
