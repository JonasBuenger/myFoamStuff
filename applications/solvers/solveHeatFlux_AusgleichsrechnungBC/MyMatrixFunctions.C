#include "MyMatrixFunctions.H"

Foam::labelList Foam::getNonZeroColumns(const scalarRectangularMatrix& A){

    labelList columns;
    bool allZero;

    for(int c=0; c<A.m(); c++){
        allZero = true;
        for(int r=0; r<A.n(); r++){
            if(A[r][c] != 0){
                allZero = false;
                break;
            }
        }
        if(!allZero)
            columns.append(c);
    }

    return columns;

}

Foam::scalarRectangularMatrix Foam::reduceColumns(const scalarRectangularMatrix& A, const labelList nonZeroColumns){

    scalarRectangularMatrix reducedA (A.n(),nonZeroColumns.size());
    int col = 0;

    forAllConstIter(labelList, nonZeroColumns, nonZeroColIter){
        for(int row=0; row<A.n(); row++){
            reducedA[row][col] = A[row][*nonZeroColIter];
        }
        col++;
    }

    return reducedA;
}

void Foam::printMatrix(const Foam::scalarRectangularMatrix& A){
    for(int row=0; row<A.n(); row++){
        for(int col=0; col<A.m(); col++){
            Info << A[row][col] << "\t";
        }
        Info << endl;
    }
}

void Foam::printMatrix(const Foam::scalarSquareMatrix& A){
    for(int row=0; row<A.n(); row++){
        for(int col=0; col<A.m(); col++){
            Info << A[row][col] << "\t";
        }
        Info << endl;
    }
}

Foam::scalarSquareMatrix Foam::convertRectToSquare(const Foam::scalarRectangularMatrix& A){

    // Gibt obere linke Quadratische Teilmatrix zurück

    int min = A.n();
    if (min > A.m())
        min = A.m();

    Foam::scalarSquareMatrix resMat(min);

    for(int row=0; row<min; row++){
        for(int col=0; col<min; col++){
            resMat[row][col] = A[row][col];
        }
    }

    return resMat;
}

Foam::scalarField Foam::convertRect2Field(const Foam::scalarRectangularMatrix& A){

    // Gibt die erste Spalte zurück
    Foam::scalarField resField(A.n());

    for(int row=0; row<A.n(); row++){
        resField[row] = A[row][0];
    }
    return resField;

}
