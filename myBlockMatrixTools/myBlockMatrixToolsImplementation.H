
#include "blockMatrixTools.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace myBlockMatrixTools
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


template<class BlockType>
void blockAdd
(
    const direction dir,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    forAll (x, i)
    {
        blockX[i](dir) += x[i];
    }
}

template<class BlockType>
void blockInsert
(
    const direction dir,
    const scalarField& x,
    Field<BlockType>& blockX
)
{
    forAll (x, i)
    {
        blockX[i](dir) = x[i];
    }
}


template<class BlockType>
void addDiag
(
    const direction dir,
    const scalarField& diag,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockB
)
{

    if (blockM.diag().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.diag().asScalar() += diag;
    }
    else if
    (
        blockM.diag().activeType() == blockCoeffBase::SCALAR
     || blockM.diag().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockDiag =
            blockM.diag().asLinear();

        forAll (diag, i)
        {
            blockDiag[i](dir) += diag[i];
        }
    }
    else if (blockM.diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<BlockType>::squareTypeField& blockDiag =
            blockM.diag().asSquare();

        forAll (diag, i)
        {
            blockDiag[i](dir, dir) += diag[i];
        }
    }

}


template<class BlockType>
void addSource
(
    const direction dir,
    const scalarField& src,
    Field<BlockType>& blockB
)
{
    blockAdd(dir, src, blockB);
}


template<class BlockType>
void addUpper
(
    const direction dir,
    const scalarField& upper,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (blockM.upper().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.upper().asScalar() += upper;
    }
    else if
    (
        blockM.upper().activeType() == blockCoeffBase::SCALAR
    || blockM.upper().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockUpper =
           blockM.upper().asLinear();
        forAll (upper, i)
        {
            blockUpper[i](dir) += upper[i];
    }
    }
    else if (blockM.upper().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<BlockType>::squareTypeField& blockUpper =
            blockM.upper().asSquare();
        forAll (upper, i)
        {
            blockUpper[i](dir, dir) += upper[i];
        }
    }
}

template<class BlockType>
void addLower
(
    const direction dir,
    const scalarField& lower,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (blockM.lower().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.lower().asScalar() += lower;
    }
    else if
    (
        blockM.lower().activeType() == blockCoeffBase::SCALAR
    || blockM.lower().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockLower =
            blockM.lower().asLinear();
        forAll (lower, i)
        {
            blockLower[i](dir) += lower[i];
        }
    }
    else if (blockM.lower().activeType() == blockCoeffBase::SQUARE)
    {
       typename CoeffField<BlockType>::squareTypeField& blockLower =
           blockM.lower().asSquare();
       forAll (lower, i)
       {
           blockLower[i](dir, dir) += lower[i];
       }
    }
}



template<class BlockType>
void insertDiag
(
    const direction dir,
    const scalarField& diag,
    BlockLduMatrix<BlockType>& blockM
)
{

    if (blockM.diag().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.diag().asScalar() = diag;
    }
    else if
    (
        blockM.diag().activeType() == blockCoeffBase::SCALAR
     || blockM.diag().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockDiag =
            blockM.diag().asLinear();

        forAll (diag, i)
        {
            blockDiag[i](dir) = diag[i];
        }
    }
    else if (blockM.diag().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<BlockType>::squareTypeField& blockDiag =
            blockM.diag().asSquare();

        forAll (diag, i)
        {
            blockDiag[i](dir, dir) = diag[i];
        }
    }

}


template<class BlockType>
void insertSource
(
    const direction dir,
    const scalarField& src,
    Field<BlockType>& blockB
)
{
    blockInsert(dir, src, blockB);
}


template<class BlockType>
void insertUpper
(
    const direction dir,
    const scalarField& upper,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (blockM.upper().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.upper().asScalar() = upper;
    }
    else if
    (
        blockM.upper().activeType() == blockCoeffBase::SCALAR
    || blockM.upper().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockUpper =
           blockM.upper().asLinear();
        forAll (upper, i)
        {
            blockUpper[i](dir) = upper[i];
    }
    }
    else if (blockM.upper().activeType() == blockCoeffBase::SQUARE)
    {
        typename CoeffField<BlockType>::squareTypeField& blockUpper =
            blockM.upper().asSquare();
        forAll (upper, i)
        {
            blockUpper[i](dir, dir) = upper[i];
        }
    }
}

template<class BlockType>
void insertLower
(
    const direction dir,
    const scalarField& lower,
    BlockLduMatrix<BlockType>& blockM
)
{
    if (blockM.lower().activeType() == blockCoeffBase::UNALLOCATED)
    {
        blockM.lower().asScalar() = lower;
    }
    else if
    (
        blockM.lower().activeType() == blockCoeffBase::SCALAR
    || blockM.lower().activeType() == blockCoeffBase::LINEAR
    )
    {
        typename CoeffField<BlockType>::linearTypeField& blockLower =
            blockM.lower().asLinear();
        forAll (lower, i)
        {
            blockLower[i](dir) = lower[i];
        }
    }
    else if (blockM.lower().activeType() == blockCoeffBase::SQUARE)
    {
       typename CoeffField<BlockType>::squareTypeField& blockLower =
           blockM.lower().asSquare();
       forAll (lower, i)
       {
           blockLower[i](dir, dir) = lower[i];
       }
    }
}

template<class BlockType>
void addCouplingSource
(
    const int dirI,
    const scalarField& source,
    Field<BlockType>& blockB
)
{
    blockAdd(dirI, source, blockB);   
}

template<class BlockType>
void addCouplingDiag
(
    const int dirI,
    const int dirJ,
    const scalarField& diag,
    BlockLduMatrix<BlockType>& blockM
)
{
    // Add off-diagonal block coefficients
    typename CoeffField<BlockType>::squareTypeField& blockDiag =
        blockM.diag().asSquare();

    // Set off-diagonal coefficient
    forAll (diag, i)
    {
        blockDiag[i](dirI, dirJ) += diag[i];
    }
}


template<class BlockType>
void addCouplingUpper
(
    const direction dirI,
    const direction dirJ,
    const scalarField& upper,
    BlockLduMatrix<BlockType>& blockM
)
{
    typename CoeffField<BlockType>::squareTypeField& blockUpper =
        blockM.upper().asSquare();

    forAll (upper, i)
    {
        blockUpper[i](dirI, dirJ) += upper[i];
    }

}


template<class BlockType>
void addCouplingLower
(
    const direction dirI,
    const direction dirJ,
    const scalarField& lower,
    BlockLduMatrix<BlockType>& blockM
)
{
    typename CoeffField<BlockType>::squareTypeField& blockLower =
        blockM.lower().asSquare();
    forAll (lower, i)
    {
        blockLower[i](dirI, dirJ) += lower[i];
    }
}

template<class BlockType>
void addVectorDiv
(
    const direction numEqn,
    const direction firstDir,
    const vectorField& diag,
    const vectorField& upper,
    const vectorField& lower,
    const scalarField& source,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockB
)
{

    myBlockMatrixTools::addCouplingDiag(numEqn,firstDir+0,diag.component(0),blockM);
    myBlockMatrixTools::addCouplingDiag(numEqn,firstDir+1,diag.component(1),blockM);
    myBlockMatrixTools::addCouplingDiag(numEqn,firstDir+2,diag.component(2),blockM);

    myBlockMatrixTools::addCouplingUpper(numEqn,firstDir+0,upper.component(0),blockM);
    myBlockMatrixTools::addCouplingUpper(numEqn,firstDir+1,upper.component(1),blockM);
    myBlockMatrixTools::addCouplingUpper(numEqn,firstDir+2,upper.component(2),blockM);

    myBlockMatrixTools::addCouplingLower(numEqn,firstDir+0,lower.component(0),blockM);
    myBlockMatrixTools::addCouplingLower(numEqn,firstDir+1,lower.component(1),blockM);
    myBlockMatrixTools::addCouplingLower(numEqn,firstDir+2,lower.component(2),blockM);

    myBlockMatrixTools::addSource(numEqn,source,blockB);

}


template<class BlockType>
void insertVectorEquation
(
    const direction firstEqn,
    fvVectorMatrix& m,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockX,
    Field<BlockType>& blockB
)
{

    if(m.hasDiag()){
        m.addBoundaryDiag(m.diag(),0);
        const scalarField& diag = m.diag();
        insertDiag  (firstEqn+0, diag,  blockM);
        insertDiag  (firstEqn+1, diag,  blockM);
        insertDiag  (firstEqn+2, diag,  blockM);
    }
    if(m.hasUpper()){
        const scalarField& upper = m.upper();
        insertUpper (firstEqn+0, upper, blockM);
        insertUpper (firstEqn+1, upper, blockM);
        insertUpper (firstEqn+2, upper, blockM);
    }
    if(m.hasLower() || ( m.symmetric() && m.hasUpper() )){
        const scalarField& lower = m.lower();
        insertLower (firstEqn+0, lower, blockM);
        insertLower (firstEqn+1, lower, blockM);
        insertLower (firstEqn+2, lower, blockM);
    }
    m.addBoundarySource(m.source(),false);
    insertSource(firstEqn+0, m.source().component(0),blockB);
    insertSource(firstEqn+1, m.source().component(1),blockB);
    insertSource(firstEqn+2, m.source().component(2),blockB);
}

template<class BlockType>
void addScalarGrad
(
    const direction firstEqn,
    const direction cmptDir,
    const vectorField& diag,
    const vectorField& upper,
    const vectorField& lower,
    const vectorField& source,
    BlockLduMatrix<BlockType>& blockM,
    Field<BlockType>& blockB
)
{

    myBlockMatrixTools::addCouplingDiag(firstEqn+0, cmptDir, diag.component(0), blockM);
    myBlockMatrixTools::addCouplingDiag(firstEqn+1, cmptDir, diag.component(1), blockM);
    myBlockMatrixTools::addCouplingDiag(firstEqn+2, cmptDir, diag.component(2), blockM);

    myBlockMatrixTools::addCouplingUpper(firstEqn+0, cmptDir, upper.component(0), blockM);
    myBlockMatrixTools::addCouplingUpper(firstEqn+1, cmptDir, upper.component(1), blockM);
    myBlockMatrixTools::addCouplingUpper(firstEqn+2, cmptDir, upper.component(2), blockM);

    myBlockMatrixTools::addCouplingLower(firstEqn+0, cmptDir, lower.component(0), blockM);
    myBlockMatrixTools::addCouplingLower(firstEqn+1, cmptDir, lower.component(1), blockM);
    myBlockMatrixTools::addCouplingLower(firstEqn+2, cmptDir, lower.component(2), blockM);

    myBlockMatrixTools::addSource(firstEqn+0, source.component(0), blockB);
    myBlockMatrixTools::addSource(firstEqn+1, source.component(1), blockB);
    myBlockMatrixTools::addSource(firstEqn+2, source.component(2), blockB);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace myBlockMatrixTools

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
