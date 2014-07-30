#include "approximations4ThetaAtBoundary.H"
#include "equalizationBasis.h"
#include "MyMatrixFunctions.H"
#include "fvCFD.H"

scalarField Foam::LeastSquareApproximation(const int type, const fvPatch& patch, const objectRegistry& db, scalarField boundaryTheta){

    // Variablen
    labelList neighbours_multiple;
    labelList neighbours_unique;
    labelList uniqueEntries;
    const fvMesh& m = patch.boundaryMesh().mesh();
    const volVectorField& c = m.C();
    const labelListList& pc = m.pointCells();
    const labelListList& cp = m.cellPoints();
    const volScalarField& internalTheta = db.lookupObject<volScalarField>("Theta");


    // LOOP OVER ALL BOUNDARY CELLS
    int face = 0;
    forAllConstIter(labelUList, patch.faceCells(), faceCell){

//        Info << "faceCell: " << *faceCell << endl;

        // BUILD LIST OF ALL NEIGHBOURING CELLS
        neighbours_multiple.clear();
        forAllConstIter(labelList,cp[*faceCell], Iter2Point){
            neighbours_multiple.append(pc[*Iter2Point]);
        }

        uniqueEntries.clear();
        Foam::uniqueOrder(neighbours_multiple, uniqueEntries);
        neighbours_unique.clear();
        forAllIter(labelList,uniqueEntries,e){
            neighbours_unique.append(neighbours_multiple[*e]);
        }

        // LOOP VARIABLES/OBJECTS
        equalizationBasis equBasis(type);
        int numBasisFcns = equBasis.dim();
        vector faceCenter = patch.Cf()[face];
        int refValues = neighbours_unique.size();
        scalarRectangularMatrix A(refValues,numBasisFcns);
        scalarRectangularMatrix b(refValues,1);

        // LOOP OVER ALL NEIGHBOURING CELLS
        int row=0;
        forAllIter(labelList,neighbours_unique, IterNeighbour){
            // Setze Koeffizienten von A und b
            for(int col=0; col<numBasisFcns; col++){
                A[row][col] = equBasis.eval(faceCenter,c[*IterNeighbour],col);
            }
            //b[row][0] = internalTheta[*IterNeighbour];
            const double f = 0.5;
            b[row][0] = f * internalTheta[*faceCell] + (1.0-f) * internalTheta[*IterNeighbour];
//            Info << "internalTheta[" << *IterNeighbour << "]: " << internalTheta[*IterNeighbour] << endl;
            row++;
        }

        // NULLSPALTEN ENTFERNEN
        labelList nonZeroColumns = getNonZeroColumns(A); // nonZeroColumns = nonZeroBasisFunctions
        A = reduceColumns(A,nonZeroColumns);
        scalarField ATb = convertRect2Field(A.T()*b);

        // GLEICHUNGSSYSTEM DEFINIEREN
        scalarSquareMatrix ATA = convertRectToSquare(A.T()*A);

        // GLEICHUNGSSYSTEM LÖSEN
        scalarField basisCoeffs(ATb.size());
        solve(basisCoeffs, ATA, ATb);

        // TEMPERATUR AM RAND SETZEN
        boundaryTheta[face] = 0;
        int indexAlpha = 0;
        forAllConstIter(labelList, nonZeroColumns, iter){
            boundaryTheta[face] += basisCoeffs[indexAlpha++]*equBasis.eval(faceCenter, faceCenter, *iter);
            //boundaryTheta[face] += basisCoeffs[indexAlpha++]*equBasis.eval(faceCenter, c[*faceCell], *iter);
        }

//        Info << "Internal: " << internalTheta[*faceCell] << endl;
//        Info << "boundary: " << boundaryTheta[face] << endl;
        face++;
    }

    return boundaryTheta;

}
