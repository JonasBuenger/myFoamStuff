#ifndef WRITECELLCENTERS_H
#define WRITECELLCENTERS_H

    volVectorField cc
    (
        IOobject
        (
            "Centers",
            runTime.constant(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero",dimless,vector(0,0,0))
    );

    cc.internalField() = mesh.C();
    cc.write();

#endif