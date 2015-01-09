/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockCoupledScalarTransportFoam

Description
    Solves two coupled transport equations in a block-coupled manner

        1) transport equation for a passive scalar
        2) diffusion only

    This resembles heat exchanging flow through a porous medium

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockLduSolvers.H"
#include "blockVectorNMatrices.H"
#include "blockMatrixTools.H"
#include "vector10/blockVector10Matrix.H"
#include "vector10.H"
#include "tensor10.H"
#include "../../../discretisation/finiteVolume/implicit/myFvm.H"
#include "../../../blockMatrix/myBlockMatrices/momentumTransportMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "writeCellCenters.H"
    #include "writeVolumes.H"

    Info<< "\nCalculating momentum transport\n" << endl;

    momentumTransportMatrix momTrans(mesh, u, p, sigma);

    bool updateOnlyRHS = false;

    momTrans.updateMatrix(updateOnlyRHS);

    updateOnlyRHS = true;

    while(runTime.loop()) //for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"
#       include "readPISOControls.H"

        momTrans.updateMatrix(updateOnlyRHS);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {

            momTrans.solve();

            momTrans.updateFields();

            momTrans.writeToFile();
        }

        runTime.write();

    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
