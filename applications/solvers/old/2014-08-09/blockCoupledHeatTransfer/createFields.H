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

    // Working coupled solution field
    Info<< "Creating field block_sT\n" << endl;
    volVector4Field block_sT
    (
        IOobject
        (
            "block_sT",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector4("zero", dimless, vector4::zero)
    );

    {

        vector4Field& blockX = block_sT.internalField();

        // insert s
        for(int i=0; i<3; i++){

            tmp<scalarField> tf = s.internalField().component(i);
            scalarField& si = tf();

            blockMatrixTools::blockInsert(i, si, blockX);
        }

        // insert Theta T
        scalarField& f = Theta.internalField();
        blockMatrixTools::blockInsert(3, f, blockX);

   }

    block_sT.write();


    label ThetaRefCell = 0;
    scalar ThetaRefValue = 0.0;
    setRefCell(Theta, mesh.solutionDict().subDict("PISO"), ThetaRefCell, ThetaRefValue);
