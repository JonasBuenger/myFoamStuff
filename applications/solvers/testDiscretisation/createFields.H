    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );


    dimensionedScalar beta 
    (
        transportProperties.lookup("beta")
    );

    dimensionedScalar corrDim
    (
        transportProperties.lookup("corrDim")
    );

    Info<< "Reading field s\n" << endl;
    volScalarField s
    (
        IOobject
        (
            "s",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field sc\n" << endl;
    volScalarField sc
    (
        IOobject
        (
            "sc",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field v\n" << endl;
    volVectorField v
    (
        IOobject
        (
            "v",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field vc\n" << endl;
    volVectorField vc
    (
        IOobject
        (
            "vc",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field sT\n" << endl;
    volSymmTensorField sT
    (
        IOobject
        (
            "sT",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field sTc\n" << endl;
    volSymmTensorField sTc
    (
        IOobject
        (
            "sTc",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Working coupled solution field
//    Info<< "Creating field X1\n" << endl;
//    volVector1Field X1
//    (
//        IOobject
//        (
//            "X1",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::NO_WRITE
//        ),
//        mesh,
//        dimensionedVector1("zero", dimless, vector1::zero)
//    );

    Info<< "Creating field X3\n" << endl;
    volVector3Field X3
    (
        IOobject
        (
            "X3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector3("zero", dimless, vector3::zero)
    );

    Info<< "Creating field X6\n" << endl;
    volVector6Field X6
    (
        IOobject
        (
            "X6",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector6("zero", dimless, vector6::zero)
    );

    surfaceScalarField phiVc
    (
        IOobject
        (
            "phiVc",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(vc) & mesh.Sf()
    );

