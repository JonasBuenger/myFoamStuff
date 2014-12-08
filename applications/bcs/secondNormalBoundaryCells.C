#include "fvCFD.H"
#include "oppositeFace.H"
#include "fvPatch.H"
#include "secondNormalBoundaryCells.H"

//template<class type>
secondNormalBoundaryCells::secondNormalBoundaryCells(const fvPatch& p):
    secondNormalFaceIsOwner(p.size()),
    secondDeltas(p.size()),
    secondFaces(p.size()),
    secondInternalCoeffs(p.size()){

    const fvMesh& m = p.boundaryMesh().mesh();
    const faceList& meshFaces = m.faces();
    const cellList cells = m.cells();
    const labelList& faceCells = p.faceCells();
    const volVectorField& centers = m.C();
    const labelList& neighbours = m.neighbour();
    const labelList& owners = m.owner();

    forAll(faceCells,i){

        const label curCellLabel = faceCells[i];
        const label curPatchFaceLabel = p.patch().start() + i;

        const label secondFace = cells[curCellLabel].opposingFaceLabel(curPatchFaceLabel,meshFaces);
        secondFaces[i] = secondFace;

        if (curCellLabel == m.owner()[secondFace]){

            secondNormalFaceIsOwner[i] = false;
            secondDeltas[i] = mag(centers[curCellLabel] - centers[neighbours[secondFace]]);

        }
        else{

            secondNormalFaceIsOwner[i] = true;
            secondDeltas[i] = mag(centers[curCellLabel] - centers[owners[secondFace]]);

        }

    }



}

//template<class type>
//const Field<type>& secondNormalBoundaryCells<type>::valueSecondInternalCoeffs(){
//    return secondInternalCoeffs;
//}

