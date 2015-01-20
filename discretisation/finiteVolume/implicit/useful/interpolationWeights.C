#include "fvCFD.H"

// upwind
// linear upwind
// linear
// limitedLinear

template<class type>
const tmp<surfaceScalarField> interpolWeights(GeometricField<type, fvPatchField, volMesh>& vf, const word& name){

    //vf.boundaryField().updateCoeffs();
    const objectRegistry& db = vf.db();

    //const volVectorField& u = db.lookupObject<volVectorField>("vc");
    const volVectorField& u = db.lookupObject<volVectorField>("u");
    //const volVectorField& u = db.lookupObject<volVectorField>("s");

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            vf.time().timeName(),
            vf.mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(u) & vf.mesh().Sf()
    );

    tmp<surfaceInterpolationScheme<type> > tinterpScheme;
    tinterpScheme = fvc::scheme<type>(phi, vf.mesh().schemesDict().interpolationScheme(name));

    const surfaceInterpolationScheme<type>& interpolationScheme = tinterpScheme();

    tmp<surfaceScalarField> tweights = interpolationScheme.weights(vf);

    //Info << "weights[15]:" << tweights()[15] << endl;

    return tweights;

}
