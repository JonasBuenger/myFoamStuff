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

#include "vector10Field.H"
#include "tensor10Field.H"
#include "fvCFD.H"

#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockLduSolvers.H"
#include "blockVectorNMatrices.H"
#include "blockMatrixTools.H"
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
typedef BlockLduPrecon<vector10> blockVector10Precon;
defineNamedTemplateTypeNameAndDebug(blockVector10Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector10Precon, dictionary);

typedef BlockNoPrecon<vector10> blockNoPreconVector10;
makeBlockPrecon(blockVector10Precon, blockNoPreconVector10);

typedef BlockDiagonalPrecon<vector10> blockDiagonalPreconVector10;
makeBlockPrecon(blockVector10Precon, blockDiagonalPreconVector10);

typedef BlockGaussSeidelPrecon<vector10> blockGaussSeidelPreconVector10;
makeBlockPrecon(blockVector10Precon, blockGaussSeidelPreconVector10);

typedef BlockCholeskyPrecon<vector10> blockCholeskyPreconVector10;
makeBlockPrecon(blockVector10Precon, blockCholeskyPreconVector10);


// Smoothers
typedef BlockLduSmoother<vector10> blockVector10Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector10Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector10Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector10> blockGaussSeidelSmootherVector10;
makeBlockSmoother(blockVector10Smoother, blockGaussSeidelSmootherVector10);

typedef BlockILUSmoother<vector10> blockILUSmootherVector10;
makeBlockSmoother(blockVector10Smoother, blockILUSmootherVector10);


// Solvers
typedef BlockLduSolver<vector10> blockVector10Solver;
defineNamedTemplateTypeNameAndDebug(blockVector10Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector10Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector10Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector10> blockDiagonalSolverVector10;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector10, 0);

typedef BlockBiCGStabSolver<vector10> blockBiCGStabSolverVector10;
makeBlockSolverTypeName(blockBiCGStabSolverVector10);
addSolverToBlockMatrix(Vector10, blockBiCGStabSolverVector10, symMatrix);
addSolverToBlockMatrix(Vector10, blockBiCGStabSolverVector10, asymMatrix);

typedef BlockCGSolver<vector10> blockCGSolverVector10;
makeBlockSolverTypeName(blockCGSolverVector10);
addSolverToBlockMatrix(Vector10, blockCGSolverVector10, symMatrix);

typedef BlockGaussSeidelSolver<vector10> blockGaussSeidelSolverVector10;
makeBlockSolverTypeName(blockGaussSeidelSolverVector10);
addSolverToBlockMatrix(Vector10, blockGaussSeidelSolverVector10, symMatrix);
addSolverToBlockMatrix(Vector10, blockGaussSeidelSolverVector10, asymMatrix);

typedef BlockGMRESSolver<vector10> blockGMRESSolverVector10;
makeBlockSolverTypeName(blockGMRESSolverVector10);
addSolverToBlockMatrix(Vector10, blockGMRESSolverVector10, symMatrix);
addSolverToBlockMatrix(Vector10, blockGMRESSolverVector10, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
