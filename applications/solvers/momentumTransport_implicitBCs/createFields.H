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
