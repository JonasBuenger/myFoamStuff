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

#include "vector3Field.H"
#include "tensor3Field.H"
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
typedef BlockLduPrecon<vector3> blockVector3Precon;
defineNamedTemplateTypeNameAndDebug(blockVector3Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector3Precon, dictionary);

typedef BlockNoPrecon<vector3> blockNoPreconVector3;
makeBlockPrecon(blockVector3Precon, blockNoPreconVector3);

typedef BlockDiagonalPrecon<vector3> blockDiagonalPreconVector3;
makeBlockPrecon(blockVector3Precon, blockDiagonalPreconVector3);

typedef BlockGaussSeidelPrecon<vector3> blockGaussSeidelPreconVector3;
makeBlockPrecon(blockVector3Precon, blockGaussSeidelPreconVector3);

typedef BlockCholeskyPrecon<vector3> blockCholeskyPreconVector3;
makeBlockPrecon(blockVector3Precon, blockCholeskyPreconVector3);


// Smoothers
typedef BlockLduSmoother<vector3> blockVector3Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector3Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector3Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector3> blockGaussSeidelSmootherVector3;
makeBlockSmoother(blockVector3Smoother, blockGaussSeidelSmootherVector3);

typedef BlockILUSmoother<vector3> blockILUSmootherVector3;
makeBlockSmoother(blockVector3Smoother, blockILUSmootherVector3);


// Solvers
typedef BlockLduSolver<vector3> blockVector3Solver;
defineNamedTemplateTypeNameAndDebug(blockVector3Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector3Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector3Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector3> blockDiagonalSolverVector3;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector3, 0);

typedef BlockBiCGStabSolver<vector3> blockBiCGStabSolverVector3;
makeBlockSolverTypeName(blockBiCGStabSolverVector3);
addSolverToBlockMatrix(Vector3, blockBiCGStabSolverVector3, symMatrix);
addSolverToBlockMatrix(Vector3, blockBiCGStabSolverVector3, asymMatrix);

typedef BlockCGSolver<vector3> blockCGSolverVector3;
makeBlockSolverTypeName(blockCGSolverVector3);
addSolverToBlockMatrix(Vector3, blockCGSolverVector3, symMatrix);

typedef BlockGaussSeidelSolver<vector3> blockGaussSeidelSolverVector3;
makeBlockSolverTypeName(blockGaussSeidelSolverVector3);
addSolverToBlockMatrix(Vector3, blockGaussSeidelSolverVector3, symMatrix);
addSolverToBlockMatrix(Vector3, blockGaussSeidelSolverVector3, asymMatrix);

typedef BlockGMRESSolver<vector3> blockGMRESSolverVector3;
makeBlockSolverTypeName(blockGMRESSolverVector3);
addSolverToBlockMatrix(Vector3, blockGMRESSolverVector3, symMatrix);
addSolverToBlockMatrix(Vector3, blockGMRESSolverVector3, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
