#ifndef WRITEVOLUMES_H
#define WRITEVOLUMES_H

    volScalarField cv
    (
        IOobject
        (
            "Volumes",
            runTime.constant(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero",dimVolume,0.0)
    );

    cv.internalField() = mesh.V();
    cv.write();


#endif 
