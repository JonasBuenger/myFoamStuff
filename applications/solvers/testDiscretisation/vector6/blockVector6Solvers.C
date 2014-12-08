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

#include "vector6Field.H"
#include "tensor6Field.H"
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
typedef BlockLduPrecon<vector6> blockVector6Precon;
defineNamedTemplateTypeNameAndDebug(blockVector6Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector6Precon, dictionary);

typedef BlockNoPrecon<vector6> blockNoPreconVector6;
makeBlockPrecon(blockVector6Precon, blockNoPreconVector6);

typedef BlockDiagonalPrecon<vector6> blockDiagonalPreconVector6;
makeBlockPrecon(blockVector6Precon, blockDiagonalPreconVector6);

typedef BlockGaussSeidelPrecon<vector6> blockGaussSeidelPreconVector6;
makeBlockPrecon(blockVector6Precon, blockGaussSeidelPreconVector6);

typedef BlockCholeskyPrecon<vector6> blockCholeskyPreconVector6;
makeBlockPrecon(blockVector6Precon, blockCholeskyPreconVector6);


// Smoothers
typedef BlockLduSmoother<vector6> blockVector6Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector6Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector6Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector6> blockGaussSeidelSmootherVector6;
makeBlockSmoother(blockVector6Smoother, blockGaussSeidelSmootherVector6);

typedef BlockILUSmoother<vector6> blockILUSmootherVector6;
makeBlockSmoother(blockVector6Smoother, blockILUSmootherVector6);


// Solvers
typedef BlockLduSolver<vector6> blockVector6Solver;
defineNamedTemplateTypeNameAndDebug(blockVector6Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector6Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector6Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector6> blockDiagonalSolverVector6;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector6, 0);

typedef BlockBiCGStabSolver<vector6> blockBiCGStabSolverVector6;
makeBlockSolverTypeName(blockBiCGStabSolverVector6);
addSolverToBlockMatrix(Vector6, blockBiCGStabSolverVector6, symMatrix);
addSolverToBlockMatrix(Vector6, blockBiCGStabSolverVector6, asymMatrix);

typedef BlockCGSolver<vector6> blockCGSolverVector6;
makeBlockSolverTypeName(blockCGSolverVector6);
addSolverToBlockMatrix(Vector6, blockCGSolverVector6, symMatrix);

typedef BlockGaussSeidelSolver<vector6> blockGaussSeidelSolverVector6;
makeBlockSolverTypeName(blockGaussSeidelSolverVector6);
addSolverToBlockMatrix(Vector6, blockGaussSeidelSolverVector6, symMatrix);
addSolverToBlockMatrix(Vector6, blockGaussSeidelSolverVector6, asymMatrix);

typedef BlockGMRESSolver<vector6> blockGMRESSolverVector6;
makeBlockSolverTypeName(blockGMRESSolverVector6);
addSolverToBlockMatrix(Vector6, blockGMRESSolverVector6, symMatrix);
addSolverToBlockMatrix(Vector6, blockGMRESSolverVector6, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
