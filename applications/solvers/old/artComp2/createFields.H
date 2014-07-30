    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    
    Info << "Reading field nu\n" << endl;
    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );
 
    Info << "Reading field beta\n" << endl;
    dimensionedScalar beta
    (
        transportProperties.lookup("beta")
    );
 

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


//#   include "createPhi.H" //brauch ich wohl nicht ???


    label ThetaRefCell = 0;
    scalar ThetaRefValue = 0.0;
    setRefCell(Theta, mesh.solutionDict().subDict("PISO"), ThetaRefCell, ThetaRefValue);
