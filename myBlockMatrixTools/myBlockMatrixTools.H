
#ifndef myBlockMatrixTools_H
#define myBlockMatrixTools_H

#include "blockLduMatrices.H"
#include "fvMatrices.H"
#include "surfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Namespace myBlockMatrixTools functions Declaration
\*---------------------------------------------------------------------------*/

namespace myBlockMatrixTools
{
    // Field operations

        //- Add field into block field
        template<class BlockType>
        void blockAdd
        (
            const direction dir,
            const scalarField& x,
            Field<BlockType>& blockX
        );

    // Diagonal block operations

        //- add matrix diagonal and source into the block system
        template<class BlockType>
        void addDiag
        (
            const direction dir,
            const scalarField& diag,
            BlockLduMatrix<BlockType>& blockM
        );

        template<class BlockType>
        void addSource
        (
            const direction dir,
            const scalarField& source,
            Field<BlockType>& blockB
        );

        // add upper and lower coefficients matrix into the block system
        template<class Type, class BlockType>
        void addUpper
        (
            const direction dir,
            const scalarField& upper,
            BlockLduMatrix<BlockType>& blockM
        );

        template<class Type, class BlockType>
        void addLower
        (
            const direction dir,
            const scalarField& lower,
            BlockLduMatrix<BlockType>& blockM
        );

    // Coupling block operations

        //- add coupling matrix diagonal and source into the block system
        template<class BlockType>
        void
        addCouplingDiag
        (
            const int dirI,
            const int dirJ,
            const scalarField& diag,
            BlockLduMatrix<BlockType>& blockM
        );

        template<class BlockType>
        void
        addCouplingSource
        (
            const int dirI,
            const scalarField& source,
            Field<BlockType>& blockB
        );

        // add coupling matrix into the block system
        template<class Type, class BlockType>
        void addCouplingUpper
        (
            const direction dirI,
            const direction dirJ,
            const scalarField& upper,
            BlockLduMatrix<BlockType>& blockM
        );

        template<class Type, class BlockType>
        void addCouplingLower
        (
            const direction dirI,
            const direction dirJ,
            const scalarField& lower,
            BlockLduMatrix<BlockType>& blockM
        );


        template<class BlockType>
        void addVectorDiv
        (
            const direction numEqn,
            const direction firstDir,
            const vectorField& diag,
            const vectorField& upper,
            const vectorField& lower,
            const vectorField& source,
            BlockLduMatrix<BlockType>& blockM,
            Field<BlockType>& blockB
        );


        template<class BlockType>
        void addScalarGrad
        (
            const direction firstEqn,
            const direction cmptDir,
            const vectorField& diag,
            const vectorField& upper,
            const vectorField& lower,
            const scalarField& source,
            BlockLduMatrix<BlockType>& blockM,
            Field<BlockType>& blockB
        );

} // End namespace myBlockMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "myBlockMatrixToolsImplementation.H"

#endif

// ************************************************************************* //
