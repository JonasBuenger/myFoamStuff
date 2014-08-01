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

#include "vector4Field.H"
#include "tensor4Field.H"
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
typedef BlockLduPrecon<vector4> blockVector4Precon;
defineNamedTemplateTypeNameAndDebug(blockVector4Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector4Precon, dictionary);

typedef BlockNoPrecon<vector4> blockNoPreconVector4;
makeBlockPrecon(blockVector4Precon, blockNoPreconVector4);

typedef BlockDiagonalPrecon<vector4> blockDiagonalPreconVector4;
makeBlockPrecon(blockVector4Precon, blockDiagonalPreconVector4);

typedef BlockGaussSeidelPrecon<vector4> blockGaussSeidelPreconVector4;
makeBlockPrecon(blockVector4Precon, blockGaussSeidelPreconVector4);

typedef BlockCholeskyPrecon<vector4> blockCholeskyPreconVector4;
makeBlockPrecon(blockVector4Precon, blockCholeskyPreconVector4);


// Smoothers
typedef BlockLduSmoother<vector4> blockVector4Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector4Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector4Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector4> blockGaussSeidelSmootherVector4;
makeBlockSmoother(blockVector4Smoother, blockGaussSeidelSmootherVector4);

typedef BlockILUSmoother<vector4> blockILUSmootherVector4;
makeBlockSmoother(blockVector4Smoother, blockILUSmootherVector4);


// Solvers
typedef BlockLduSolver<vector4> blockVector4Solver;
defineNamedTemplateTypeNameAndDebug(blockVector4Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector4Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector4Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector4> blockDiagonalSolverVector4;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector4, 0);

typedef BlockBiCGStabSolver<vector4> blockBiCGStabSolverVector4;
makeBlockSolverTypeName(blockBiCGStabSolverVector4);
addSolverToBlockMatrix(Vector4, blockBiCGStabSolverVector4, symMatrix);
addSolverToBlockMatrix(Vector4, blockBiCGStabSolverVector4, asymMatrix);

typedef BlockCGSolver<vector4> blockCGSolverVector4;
makeBlockSolverTypeName(blockCGSolverVector4);
addSolverToBlockMatrix(Vector4, blockCGSolverVector4, symMatrix);

typedef BlockGaussSeidelSolver<vector4> blockGaussSeidelSolverVector4;
makeBlockSolverTypeName(blockGaussSeidelSolverVector4);
addSolverToBlockMatrix(Vector4, blockGaussSeidelSolverVector4, symMatrix);
addSolverToBlockMatrix(Vector4, blockGaussSeidelSolverVector4, asymMatrix);

typedef BlockGMRESSolver<vector4> blockGMRESSolverVector4;
makeBlockSolverTypeName(blockGMRESSolverVector4);
addSolverToBlockMatrix(Vector4, blockGMRESSolverVector4, symMatrix);
addSolverToBlockMatrix(Vector4, blockGMRESSolverVector4, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
