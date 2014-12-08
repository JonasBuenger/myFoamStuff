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
#include "vector9/blockVector9Matrix.H"
#include "vector9.H"
#include "tensor9.H"
#include "../../../discretisation/finiteVolume/myFvm.H"
#include "../../../blockMatrix/myBlockMatrices/momentumTransportTraceFreeMatrix.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating momentum transport\n" << endl;

    momentumTransportTraceFreeMatrix momTrans(mesh, u, p, sigma);

    while(runTime.loop()) //for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {

            momTrans.updateMatrix();

            //momTrans.displayBlock(1);

            //momTrans.displayCurrVal(56);

            /*
            tmp<volTensorField> tgradU = fvc::grad(u);
            volTensorField gradU = tgradU();

            Info << gradU.internalField()[56] << endl;

            tmp<volVectorField> tlaplacianSigma = fvc::div(sigma);
            volVectorField& laplacianSigma = tlaplacianSigma();

            Info << laplacianSigma.internalField()[56] << endl;

            tmp<volTensorField> tgradU = fvc::grad(u);
            volTensorField gradU = tgradU();

            Info << gradU.internalField()[56] << endl;

            tmp<volVectorField> tgradP = fvc::grad(p);
            volVectorField& gradP = tgradP();

            Info << gradP.internalField()[56] << endl;

            tmp<volScalarField> tdivU = fvc::div(u);
            volScalarField divU = tdivU();

            Info << divU.internalField()[56] << endl;
            */


            momTrans.solve();

            momTrans.updateFields();

        }

        runTime.write();
        Info << "tr(sigma[15]): " << tr(sigma[15]) << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
