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

    Info<< "Reading boundaryConstants\n" << endl;

    IOdictionary boundaryConstants
    (
        IOobject
        (
            "boundaryConstants",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

/*
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
*/

    Info<< "Reading field Theta\n" << endl;
    volScalarField Theta 
    (
        IOobject
        (
            "Theta",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field s\n" << endl;
    volVectorField s
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

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field u\n" << endl;
    volVectorField u
    (
        IOobject
        (
            "u",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field sigma\n" << endl;
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
