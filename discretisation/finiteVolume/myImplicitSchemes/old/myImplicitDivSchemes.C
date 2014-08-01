#include "myImplicitDivSchemes.H"


tmp<fvScalarMatrix>
Foam::myImplicitDivSchemes::coefficients(const volVectorField& vf, int numComponent)
{

    Info << "Foam::myImplicitDivSchemes::coefficients" << endl;

    const fvMesh& mesh = vf.mesh();

    tmp<fvScalarMatrix > tfvm
    (
        new fvScalarMatrix
        (
            vf.component(numComponent)(),
            vf.dimensions()
        )
    );

    fvScalarMatrix& fvm = tfvm();

    //fvm.psi() = vf;

    // INTERPOLATION SCHEME
    tmp<surfaceInterpolationScheme<scalar> > tinterpScheme_(
            surfaceInterpolationScheme<scalar>::New(
                    vf.mesh(),
                    mesh.schemesDict().interpolationScheme(vf.name())
            )
    );

    tmp<surfaceScalarField> tweights = tinterpScheme_().weights(vf.component(numComponent));
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

        diag[o] +=  s*w;
        diag[n] -=  s*(1-w);
        lower[i] = -s*w;
        upper[i] =  s*(1-w);

    }

    // CONTRIBUTION BOUNDARY FIELD ...      (update coeffs in solver!!!)
    //      vfpatchFace = ic * vfboundaryCell  +  bc
    // vf.boundaryField().updateCoeffs();

    //fvm.source() = vf.oldTime().internalField().component(numComponent);

    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
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

}
