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

    //Info<< "Creating field block_sT\n" << endl;

    /*
    volVector10Field blockX
    (
       IOobject
       (
           "blockX",
           runTime.timeName(),
           mesh,
           IOobject::NO_READ,
           IOobject::NO_WRITE
       ),
       mesh,
       dimensionedVector10("zero", dimless, vector10::zero)
    );
    */
    //label ThetaRefCell = 0;
    //scalar ThetaRefValue = 0.0;
    //setRefCell(Theta, mesh.solutionDict().subDict("PISO"), ThetaRefCell, ThetaRefValue);
