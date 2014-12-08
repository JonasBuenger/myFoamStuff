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

#include "vector14Field.H"
#include "tensor14Field.H"
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
typedef BlockLduPrecon<vector14> blockVector14Precon;
defineNamedTemplateTypeNameAndDebug(blockVector14Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector14Precon, dictionary);

typedef BlockNoPrecon<vector14> blockNoPreconVector14;
makeBlockPrecon(blockVector14Precon, blockNoPreconVector14);

typedef BlockDiagonalPrecon<vector14> blockDiagonalPreconVector14;
makeBlockPrecon(blockVector14Precon, blockDiagonalPreconVector14);

typedef BlockGaussSeidelPrecon<vector14> blockGaussSeidelPreconVector14;
makeBlockPrecon(blockVector14Precon, blockGaussSeidelPreconVector14);

typedef BlockCholeskyPrecon<vector14> blockCholeskyPreconVector14;
makeBlockPrecon(blockVector14Precon, blockCholeskyPreconVector14);


// Smoothers
typedef BlockLduSmoother<vector14> blockVector14Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector14Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector14Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector14> blockGaussSeidelSmootherVector14;
makeBlockSmoother(blockVector14Smoother, blockGaussSeidelSmootherVector14);

typedef BlockILUSmoother<vector14> blockILUSmootherVector14;
makeBlockSmoother(blockVector14Smoother, blockILUSmootherVector14);


// Solvers
typedef BlockLduSolver<vector14> blockVector14Solver;
defineNamedTemplateTypeNameAndDebug(blockVector14Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector14Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector14Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector14> blockDiagonalSolverVector14;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector14, 0);

typedef BlockBiCGStabSolver<vector14> blockBiCGStabSolverVector14;
makeBlockSolverTypeName(blockBiCGStabSolverVector14);
addSolverToBlockMatrix(Vector14, blockBiCGStabSolverVector14, symMatrix);
addSolverToBlockMatrix(Vector14, blockBiCGStabSolverVector14, asymMatrix);

typedef BlockCGSolver<vector14> blockCGSolverVector14;
makeBlockSolverTypeName(blockCGSolverVector14);
addSolverToBlockMatrix(Vector14, blockCGSolverVector14, symMatrix);

typedef BlockGaussSeidelSolver<vector14> blockGaussSeidelSolverVector14;
makeBlockSolverTypeName(blockGaussSeidelSolverVector14);
addSolverToBlockMatrix(Vector14, blockGaussSeidelSolverVector14, symMatrix);
addSolverToBlockMatrix(Vector14, blockGaussSeidelSolverVector14, asymMatrix);

typedef BlockGMRESSolver<vector14> blockGMRESSolverVector14;
makeBlockSolverTypeName(blockGMRESSolverVector14);
addSolverToBlockMatrix(Vector14, blockGMRESSolverVector14, symMatrix);
addSolverToBlockMatrix(Vector14, blockGMRESSolverVector14, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
