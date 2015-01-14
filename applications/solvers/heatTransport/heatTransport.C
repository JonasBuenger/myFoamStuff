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
#include "vector4/blockVector4Matrix.H"
#include "vector4.H"
#include "tensor4.H"
#include "../../../discretisation/finiteVolume/implicit/myFvm.H"
//#include "../../../discretisation/finiteVolume/myFvm.H"
#include "../../../blockMatrix/myBlockMatrices/heatTransportMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating momentum transport\n" << endl;

    #include "writeCellCenters.H"
    #include "writeVolumes.H"

    heatTransportMatrix heatTrans(mesh, s, Theta);

    while(runTime.loop()) //for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {

/*
            // Update s
            fvVectorMatrix sEqn
            (
                 fvm::ddt(s)
               - fvm::laplacian(nu/2,s)
               + corrDim*s
            );

            solve(sEqn == -fvc::grad(Theta));

            // Update Theta
            fvScalarMatrix ThetaEqn
            (
                 fvm::ddt(Theta)
            );
            solve(ThetaEqn == -fvc::div(s));

            // Output
            runTime.write();

            Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
                 << "  ClockTime = " << runTime.elapsedClockTime() << " s"
                 << nl << endl;
*/

            heatTrans.updateMatrix();

            heatTrans.solve();

            heatTrans.updateFields();

        }

        runTime.write();

    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
