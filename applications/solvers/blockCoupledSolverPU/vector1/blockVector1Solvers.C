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

\*---------------------------------------------------------------------------*/

#include "vector1Field.H"
#include "tensor1Field.H"
#include "ExpandTensorN.H"
#include "ExpandTensorNField.H"
#include "blockLduMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "blockLduPrecons.H"
#include "BlockNoPrecon.H"
#include "blockDiagonalPrecons.H"
#include "blockGaussSeidelPrecons.H"
#include "BlockCholeskyPrecon.H"

#include "blockLduSmoothers.H"
#include "blockGaussSeidelSmoothers.H"
#include "BlockILUSmoother.H"

#include "blockLduSolvers.H"
#include "BlockDiagonalSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"
#include "BlockGaussSeidelSolver.H"
#include "BlockGMRESSolver.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Preconditioners
typedef BlockLduPrecon<vector1> blockVector1Precon;
defineNamedTemplateTypeNameAndDebug(blockVector1Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector1Precon, dictionary);

typedef BlockNoPrecon<vector1> blockNoPreconVector1;
makeBlockPrecon(blockVector1Precon, blockNoPreconVector1);

typedef BlockDiagonalPrecon<vector1> blockDiagonalPreconVector1;
makeBlockPrecon(blockVector1Precon, blockDiagonalPreconVector1);

typedef BlockGaussSeidelPrecon<vector1> blockGaussSeidelPreconVector1;
makeBlockPrecon(blockVector1Precon, blockGaussSeidelPreconVector1);

typedef BlockCholeskyPrecon<vector1> blockCholeskyPreconVector1;
makeBlockPrecon(blockVector1Precon, blockCholeskyPreconVector1);


// Smoothers
typedef BlockLduSmoother<vector1> blockVector1Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector1Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector1Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector1> blockGaussSeidelSmootherVector1;
makeBlockSmoother(blockVector1Smoother, blockGaussSeidelSmootherVector1);

typedef BlockILUSmoother<vector1> blockILUSmootherVector1;
makeBlockSmoother(blockVector1Smoother, blockILUSmootherVector1);


// Solvers
typedef BlockLduSolver<vector1> blockVector1Solver;
defineNamedTemplateTypeNameAndDebug(blockVector1Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector1Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector1Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector1> blockDiagonalSolverVector1;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector1, 0);

typedef BlockBiCGStabSolver<vector1> blockBiCGStabSolverVector1;
makeBlockSolverTypeName(blockBiCGStabSolverVector1);
addSolverToBlockMatrix(Vector1, blockBiCGStabSolverVector1, symMatrix);
addSolverToBlockMatrix(Vector1, blockBiCGStabSolverVector1, asymMatrix);

typedef BlockCGSolver<vector1> blockCGSolverVector1;
makeBlockSolverTypeName(blockCGSolverVector1);
addSolverToBlockMatrix(Vector1, blockCGSolverVector1, symMatrix);

typedef BlockGaussSeidelSolver<vector1> blockGaussSeidelSolverVector1;
makeBlockSolverTypeName(blockGaussSeidelSolverVector1);
addSolverToBlockMatrix(Vector1, blockGaussSeidelSolverVector1, symMatrix);
addSolverToBlockMatrix(Vector1, blockGaussSeidelSolverVector1, asymMatrix);

typedef BlockGMRESSolver<vector1> blockGMRESSolverVector1;
makeBlockSolverTypeName(blockGMRESSolverVector1);
addSolverToBlockMatrix(Vector1, blockGMRESSolverVector1, symMatrix);
addSolverToBlockMatrix(Vector1, blockGMRESSolverVector1, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
