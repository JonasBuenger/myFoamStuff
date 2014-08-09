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

#include "myBlockMatrixTools.H"

#include "myImplicitSchemes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating heat transport\n" << endl;

    while(runTime.loop()) //for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {

            /**************************************************************
             * CONSTRUCT BLOCK-MATRIX M AND BLOCK-SOURCE-VECTOR S
             *************************************************************/

            // block-structure
            //
            // unknown vector:  s.x()
            //                  s.y()
            //                  s.z()
            //                  Theta
            //
            // matrix:          0- heatTransferEqn.x()
            //                  1- heatTransferEqn.y()
            //                  2- heatTransferEqn.z()
            //                  3- continuityEqn

            // ... block matrix blockM
            BlockLduMatrix<vector4> blockM(mesh);

            blockM.diag().asSquare() = tensor4::zero;
            blockM.upper().asSquare() = tensor4::zero;
            blockM.lower().asSquare() = tensor4::zero;

            // ... block source vector blockB
            Field<vector4> blockB(mesh.nCells(), vector4::zero);
            blockB = vector4::zero;


            /**************************************************************
             * SET MATRIX COEFFICIENTS
             *************************************************************/

            // HEAT TRANSFER ...

            // generate first part using standard FOAM-ext

            fvVectorMatrix sEqn
            (
                  fvm::ddt(s)
                - fvm::laplacian(nu,s)
                + fvm::Sp(corrDim,s)
            //  + fvc::grad(Theta)
            //  == 0
            );

            myBlockMatrixTools::insertVectorEquation(0, sEqn, blockM, block_sT, blockB);

            myImplicitGradScheme ThetaGrad(Theta);

            myBlockMatrixTools::addScalarGrad(0, 3,
                                              ThetaGrad.coeffsDiag(),
                                              ThetaGrad.coeffsUpper(),
                                              ThetaGrad.coeffsLower(),
                                              ThetaGrad.coeffsSource(),
                                              blockM,
                                              blockB);

            // CONTINUITY EQUATION ...

            // generate first part using standard FOAM-ext
            //Theta.boundaryField().updateCoeffs();

            tmp<volScalarField> tA = sEqn.A();
            volScalarField& A = tA();

            tmp<volVectorField> texp = fvc::grad(Theta);
            volVectorField& exp = texp();
            tmp<volVectorField> texp2 = exp/A;
            volVectorField exp2 = texp2();

            fvScalarMatrix ThetaEqn
            (
                  (1/beta) *
                  fvm::ddt(Theta)
            //    - fvm::laplacian(1/A, Theta)
            //      ==
            //    - fvc::div(exp2)
            //  + fvc::div(s)
            //    == 0
            );

            blockMatrixTools::insertEquation(3, ThetaEqn, blockM, block_sT, blockB);

            myImplicitDivScheme sDiv(s);

            myBlockMatrixTools::addVectorDiv(3, 0,
                                             sDiv.coeffsDiag(),
                                             sDiv.coeffsUpper(),
                                             sDiv.coeffsLower(),
                                             sDiv.coeffsSource(),
                                             blockM, blockB);

            /**************************************************************
             * SOLVE SYSTEM
             *************************************************************/

            BlockSolverPerformance<vector4> solverPerf =
                BlockLduSolver<vector4>::New
                (
                    block_sT.name(),
                    blockM,
                    mesh.solutionDict().solver(block_sT.name())
                )->solve(block_sT, blockB);

            solverPerf.print();


            /**************************************************************
             * RETRIEVE SOLUTION
             *************************************************************/

            scalarField sx(s.internalField().size(),pTraits<scalar>::zero);
            scalarField sy(s.internalField().size(),pTraits<scalar>::zero);
            scalarField sz(s.internalField().size(),pTraits<scalar>::zero);

            blockMatrixTools::blockRetrieve(0, sx, block_sT);
            blockMatrixTools::blockRetrieve(1, sy, block_sT);
            blockMatrixTools::blockRetrieve(2, sz, block_sT);

            s.internalField().replace(0,sx);
            s.internalField().replace(1,sy);
            s.internalField().replace(2,sz);

            blockMatrixTools::blockRetrieve(3, Theta.internalField(), block_sT);

            // Theta.relax();

            s.correctBoundaryConditions();              // ruft evaluate() --> updated_ = false auf
            Theta.correctBoundaryConditions();;         // ruft evaluate() --> updated_ = false auf

        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
