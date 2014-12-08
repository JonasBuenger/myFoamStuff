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

#include "vector9Field.H"
#include "tensor9Field.H"
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
typedef BlockLduPrecon<vector9> blockVector9Precon;
defineNamedTemplateTypeNameAndDebug(blockVector9Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector9Precon, dictionary);

typedef BlockNoPrecon<vector9> blockNoPreconVector9;
makeBlockPrecon(blockVector9Precon, blockNoPreconVector9);

typedef BlockDiagonalPrecon<vector9> blockDiagonalPreconVector9;
makeBlockPrecon(blockVector9Precon, blockDiagonalPreconVector9);

typedef BlockGaussSeidelPrecon<vector9> blockGaussSeidelPreconVector9;
makeBlockPrecon(blockVector9Precon, blockGaussSeidelPreconVector9);

typedef BlockCholeskyPrecon<vector9> blockCholeskyPreconVector9;
makeBlockPrecon(blockVector9Precon, blockCholeskyPreconVector9);


// Smoothers
typedef BlockLduSmoother<vector9> blockVector9Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector9Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector9Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector9> blockGaussSeidelSmootherVector9;
makeBlockSmoother(blockVector9Smoother, blockGaussSeidelSmootherVector9);

typedef BlockILUSmoother<vector9> blockILUSmootherVector9;
makeBlockSmoother(blockVector9Smoother, blockILUSmootherVector9);


// Solvers
typedef BlockLduSolver<vector9> blockVector9Solver;
defineNamedTemplateTypeNameAndDebug(blockVector9Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector9Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector9Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector9> blockDiagonalSolverVector9;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector9, 0);

typedef BlockBiCGStabSolver<vector9> blockBiCGStabSolverVector9;
makeBlockSolverTypeName(blockBiCGStabSolverVector9);
addSolverToBlockMatrix(Vector9, blockBiCGStabSolverVector9, symMatrix);
addSolverToBlockMatrix(Vector9, blockBiCGStabSolverVector9, asymMatrix);

typedef BlockCGSolver<vector9> blockCGSolverVector9;
makeBlockSolverTypeName(blockCGSolverVector9);
addSolverToBlockMatrix(Vector9, blockCGSolverVector9, symMatrix);

typedef BlockGaussSeidelSolver<vector9> blockGaussSeidelSolverVector9;
makeBlockSolverTypeName(blockGaussSeidelSolverVector9);
addSolverToBlockMatrix(Vector9, blockGaussSeidelSolverVector9, symMatrix);
addSolverToBlockMatrix(Vector9, blockGaussSeidelSolverVector9, asymMatrix);

typedef BlockGMRESSolver<vector9> blockGMRESSolverVector9;
makeBlockSolverTypeName(blockGMRESSolverVector9);
addSolverToBlockMatrix(Vector9, blockGMRESSolverVector9, symMatrix);
addSolverToBlockMatrix(Vector9, blockGMRESSolverVector9, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
