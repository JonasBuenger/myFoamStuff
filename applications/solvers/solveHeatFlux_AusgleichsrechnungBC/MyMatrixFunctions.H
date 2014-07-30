#ifndef MYMATRIXFUNCTIONS_H
#define MYMATRIXFUNCTIONS_H

#include "simpleMatrix.H"

namespace Foam{

labelList getNonZeroColumns(const scalarRectangularMatrix& A);
scalarRectangularMatrix reduceColumns(const scalarRectangularMatrix& A, const labelList nonZeroColumns);
void printMatrix(const scalarRectangularMatrix& A);
void printMatrix(const scalarSquareMatrix& A);
scalarSquareMatrix convertRectToSquare(const scalarRectangularMatrix& A);
scalarField convertRect2Field(const scalarRectangularMatrix& A);

}

#endif // MYMATRIXFUNCTIONS_H
