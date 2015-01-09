#include"fvCFD.H"
#include "../useful/interpolationWeights.H"

namespace myFvm{

template<class vectorType>
void gaussGrad(volScalarField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    gaussGrad(vf,1,M,B,dirI,dirJ);

//    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

//    const word name("grad(" + vf.name() + ")");
//    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
//    const surfaceScalarField& weights = tweights();

//    // contribution of internal Field ...
//    for(int i=0; i<vf.mesh().owner().size(); i++)
//    {
//        int o = vf.mesh().owner()[i];
//        int n = vf.mesh().neighbour()[i];
//        scalar w = weights.internalField()[i];
//        vector s = vf.mesh().Sf()[i];

//        diag[o](dirI + 0, dirJ) += w * s.x();
//        diag[o](dirI + 1, dirJ) += w * s.y();
//        diag[o](dirI + 2, dirJ) += w * s.z();

//        diag[n](dirI + 0, dirJ) -= (1-w) * s.x();
//        diag[n](dirI + 1, dirJ) -= (1-w) * s.y();
//        diag[n](dirI + 2, dirJ) -= (1-w) * s.z();

//        upper[i](dirI + 0, dirJ) += (1-w) * s.x();
//        upper[i](dirI + 1, dirJ) += (1-w) * s.y();
//        upper[i](dirI + 2, dirJ) += (1-w) * s.z();

//        lower[i](dirI + 0, dirJ) -= w * s.x();
//        lower[i](dirI + 1, dirJ) -= w * s.y();
//        lower[i](dirI + 2, dirJ) -= w * s.z();

//    }

//    // CONTRIBUTION BOUNDARY FIELD ...
//    //
//    //      vfboundFace = ic * vfboundCell  +  bc

//    vf.boundaryField().updateCoeffs();
//    forAll(vf.boundaryField(), patchI){ // loop all patches

//        const fvPatchField<scalar>& patchVolField = vf.boundaryField()[patchI];
//        const fvPatch& patch = patchVolField.patch();   // reference to patch
//        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

//        tmp<vectorField> tsn = patch.Sf();
//        const vectorField sn = tsn();   // patch normals

//        if (patchVolField.type() == "implicitExtrapolation"){

//            tmp<scalarField > ticd = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//            tmp<scalarField > ticu = patchVolField.valueInternalCoeffsUpper(weightsPatchVolField);
//            tmp<scalarField > ticl = patchVolField.valueInternalCoeffsLower(weightsPatchVolField);
//            const scalarField& icd = ticd();      // internal coefficient
//            const scalarField& icu = ticu();
//            const scalarField& icl = ticl();
//            const List<label>& secondFaces = patchVolField.secondFaces();

//            forAll(patchVolField, faceI){       // loop all faces

//                label c = patch.faceCells()[faceI];     // boundary cell
//                const label secondFace = secondFaces[faceI];

//                diag[c](dirI + 0, dirJ) += icd[faceI]*sn[faceI].x();
//                diag[c](dirI + 1, dirJ) += icd[faceI]*sn[faceI].y();
//                diag[c](dirI + 2, dirJ) += icd[faceI]*sn[faceI].z();

//                upper[secondFace](dirI + 0, dirJ) += icu[faceI]*sn[faceI].x();
//                upper[secondFace](dirI + 1, dirJ) += icu[faceI]*sn[faceI].y();
//                upper[secondFace](dirI + 2, dirJ) += icu[faceI]*sn[faceI].z();

//                lower[secondFace](dirI + 0, dirJ) += icl[faceI]*sn[faceI].x();
//                lower[secondFace](dirI + 1, dirJ) += icl[faceI]*sn[faceI].y();
//                lower[secondFace](dirI + 2, dirJ) += icl[faceI]*sn[faceI].z();

//            }

//        }
//        else{

//            tmp<scalarField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//            tmp<scalarField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
//            const scalarField& ic = tic();     // internal coefficient
//            const scalarField& bc = tbc();     // boundary coefficient

//            forAll(patchVolField, faceI){   // loop all faces

//                label c = patch.faceCells()[faceI];     // boundary cell

//                diag[c](dirI + 0, dirJ) += ic[faceI]*sn[faceI].x();
//                diag[c](dirI + 1, dirJ) += ic[faceI]*sn[faceI].y();
//                diag[c](dirI + 2, dirJ) += ic[faceI]*sn[faceI].z();

//                B[c](dirI + 0) -= bc[faceI]*sn[faceI].x();
//                B[c](dirI + 1) -= bc[faceI]*sn[faceI].y();
//                B[c](dirI + 2) -= bc[faceI]*sn[faceI].z();

//            }
//        }


//    }

}

template<class vectorType>
void gaussSymmGrad(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    gaussSymmGrad(vf,1,M,B,dirI,dirJ);

//    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

//    const word name("grad(" + vf.name() + ")");
//    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
//    const surfaceScalarField& weights = tweights();

//    const int X = 0;
//    const int Y = 1;
//    const int Z = 2;

//    // contribution of internal Field ...
//    for(int i=0; i<vf.mesh().owner().size(); i++)
//    {
//        int o = vf.mesh().owner()[i];
//        int n = vf.mesh().neighbour()[i];
//        scalar w = weights.internalField()[i];
//        vector s = vf.mesh().Sf()[i];

//        // xx
//        diag[o](dirI + 0, dirJ + X) += s.x() * w;
//        // xy
//        diag[o](dirI + 1, dirJ + X) += 0.5 * s.y() * w;
//        diag[o](dirI + 1, dirJ + Y) += 0.5 * s.x() * w;
//        // xz
//        diag[o](dirI + 2, dirJ + X) += 0.5 * s.z() * w;
//        diag[o](dirI + 2, dirJ + Z) += 0.5 * s.x() * w;
//        // yy
//        diag[o](dirI + 3, dirJ + Y) += s.y() * w;
//        // yz
//        diag[o](dirI + 4, dirJ + Y) += 0.5 * s.z() * w;
//        diag[o](dirI + 4, dirJ + Z) += 0.5 * s.y() * w;
//        // zz
//        diag[o](dirI + 5, dirJ + Z) += s.z() * w;

//        // xx
//        diag[n](dirI + 0, dirJ + X) -= s.x() * (1-w);
//        // xy
//        diag[n](dirI + 1, dirJ + X) -= 0.5 * s.y() * (1-w);
//        diag[n](dirI + 1, dirJ + Y) -= 0.5 * s.x() * (1-w);
//        // xz
//        diag[n](dirI + 2, dirJ + X) -= 0.5 * s.z() * (1-w);
//        diag[n](dirI + 2, dirJ + Z) -= 0.5 * s.x() * (1-w);
//        // yy
//        diag[n](dirI + 3, dirJ + Y) -= s.y() * (1-w);
//        // yz
//        diag[n](dirI + 4, dirJ + Y) -= 0.5 * s.z() * (1-w);
//        diag[n](dirI + 4, dirJ + Z) -= 0.5 * s.y() * (1-w);
//        // zz
//        diag[n](dirI + 5, dirJ + Z) -= s.z() * (1-w);

//        // xx
//        upper[i](dirI + 0, dirJ + X) += s.x() * (1-w);
//        // xy
//        upper[i](dirI + 1, dirJ + X) += 0.5 * s.y() * (1-w);
//        upper[i](dirI + 1, dirJ + Y) += 0.5 * s.x() * (1-w);
//        // xz
//        upper[i](dirI + 2, dirJ + X) += 0.5 * s.z() * (1-w);
//        upper[i](dirI + 2, dirJ + Z) += 0.5 * s.x() * (1-w);
//        // yy
//        upper[i](dirI + 3, dirJ + Y) += s.y() * (1-w);
//        // yz
//        upper[i](dirI + 4, dirJ + Y) += 0.5 * s.z() * (1-w);
//        upper[i](dirI + 4, dirJ + Z) += 0.5 * s.y() * (1-w);
//        // zz
//        upper[i](dirI + 5, dirJ + Z) += s.z() * (1-w);

//        // xx
//        lower[i](dirI + 0, dirJ + X) -= s.x() * w;
//        // xy
//        lower[i](dirI + 1, dirJ + X) -= 0.5 * s.y() * w;
//        lower[i](dirI + 1, dirJ + Y) -= 0.5 * s.x() * w;
//        // xz
//        lower[i](dirI + 2, dirJ + X) -= 0.5 * s.z() * w;
//        lower[i](dirI + 2, dirJ + Z) -= 0.5 * s.x() * w;
//        // yy
//        lower[i](dirI + 3, dirJ + Y) -= s.y() * w;
//        // yz
//        lower[i](dirI + 4, dirJ + Y) -= 0.5 * s.z() * w;
//        lower[i](dirI + 4, dirJ + Z) -= 0.5 * s.y() * w;
//        // zz
//        lower[i](dirI + 5, dirJ + Z) -= s.z() * w;

//    }

//    // CONTRIBUTION BOUNDARY FIELD ...
//    //
//    //      vfboundFace = ic * vfboundCell  +  bc

//    vf.boundaryField().updateCoeffs();

//    forAll(vf.boundaryField(), patchI){ // loop all patches

//        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
//        const fvPatch& patch = patchVolField.patch();   // reference to patch
//        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

//        tmp<vectorField> tsn = patch.Sf();
//        const vectorField sn = tsn();   // patch normals

//        if (patchVolField.type() == "implicitExtrapolation"){

//            forAll(patchVolField, faceI){       // loop all faces

//                tmp<vectorField > ticd = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//                tmp<vectorField > ticu = patchVolField.valueInternalCoeffsUpper(weightsPatchVolField);
//                tmp<vectorField > ticl = patchVolField.valueInternalCoeffsLower(weightsPatchVolField);
//                const vectorField& icd = ticd();      // internal coefficient
//                const vectorField& icu = ticu();
//                const vectorField& icl = ticl();
//                const List<label>& secondFaces = patchVolField.secondFaces();

//                label c = patch.faceCells()[faceI];     // boundary cell
//                const label secondFace = secondFaces[faceI];
//                const vector s = sn[faceI];

//                diag[c](dirI + 0, dirJ + X) += icd[faceI].x() * s.x();      // xx
//                diag[c](dirI + 1, dirJ + X) += 0.5 * icd[faceI].x() * s.y();// xy
//                diag[c](dirI + 1, dirJ + Y) += 0.5 * icd[faceI].y() * s.x();// xy
//                diag[c](dirI + 2, dirJ + X) += 0.5 * icd[faceI].x() * s.z();// xz
//                diag[c](dirI + 2, dirJ + Z) += 0.5 * icd[faceI].z() * s.x();// xz
//                diag[c](dirI + 3, dirJ + Y) += icd[faceI].y() * s.y();      // yy
//                diag[c](dirI + 4, dirJ + Y) += 0.5 * icd[faceI].y() * s.z();// yz
//                diag[c](dirI + 4, dirJ + Z) += 0.5 * icd[faceI].z() * s.y();// yz
//                diag[c](dirI + 5, dirJ + Z) += icd[faceI].z() * s.z();      // zz

//                if(icu[faceI] != vector(0,0,0)){
//                    upper[secondFace](dirI + 0, dirJ + X) += icu[faceI].x() * s.x();      // xx
//                    upper[secondFace](dirI + 1, dirJ + X) += 0.5 * icu[faceI].x() * s.y();// xy
//                    upper[secondFace](dirI + 1, dirJ + Y) += 0.5 * icu[faceI].y() * s.x();// xy
//                    upper[secondFace](dirI + 2, dirJ + X) += 0.5 * icu[faceI].x() * s.z();// xz
//                    upper[secondFace](dirI + 2, dirJ + Z) += 0.5 * icu[faceI].z() * s.x();// xz
//                    upper[secondFace](dirI + 3, dirJ + Y) += icu[faceI].y() * s.y();      // yy
//                    upper[secondFace](dirI + 4, dirJ + Y) += 0.5 * icu[faceI].y() * s.z();// yz
//                    upper[secondFace](dirI + 4, dirJ + Z) += 0.5 * icu[faceI].z() * s.y();// yz
//                    upper[secondFace](dirI + 5, dirJ + Z) += icu[faceI].z() * s.z();      // zz
//                }
//                if(icl[faceI] != vector(0,0,0)){
//                    lower[secondFace](dirI + 0, dirJ + X) += icl[faceI].x() * s.x();      // xx
//                    lower[secondFace](dirI + 1, dirJ + X) += 0.5 * icl[faceI].x() * s.y();// xy
//                    lower[secondFace](dirI + 1, dirJ + Y) += 0.5 * icl[faceI].y() * s.x();// xy
//                    lower[secondFace](dirI + 2, dirJ + X) += 0.5 * icl[faceI].x() * s.z();// xz
//                    lower[secondFace](dirI + 2, dirJ + Z) += 0.5 * icl[faceI].z() * s.x();// xz
//                    lower[secondFace](dirI + 3, dirJ + Y) += icl[faceI].y() * s.y();      // yy
//                    lower[secondFace](dirI + 4, dirJ + Y) += 0.5 * icl[faceI].y() * s.z();// yz
//                    lower[secondFace](dirI + 4, dirJ + Z) += 0.5 * icl[faceI].z() * s.y();// yz
//                    lower[secondFace](dirI + 5, dirJ + Z) += icl[faceI].z() * s.z();      // zz
//                }

//            }

//        }
//        else if (patchVolField.type() == "myImplicitVelocity"){

//            tmp<tensorField > ticds = patchVolField.valueInternalCoeffsTensor(weightsPatchVolField);
//            tmp<tensorField > ticus = patchVolField.valueInternalCoeffsUpperTensor(weightsPatchVolField);
//            tmp<tensorField > ticls = patchVolField.valueInternalCoeffsLowerTensor(weightsPatchVolField);
//            const tensorField& icds = ticds();      // internal coefficient
//            const tensorField& icus = ticus();
//            const tensorField& icls = ticls();
//            const List<label>& secondFaces = patchVolField.secondFaces();

//            forAll(patchVolField, faceI){       // loop all faces

//                label c = patch.faceCells()[faceI];     // boundary cell
//                const label secondFace = secondFaces[faceI];
//                const vector s = sn[faceI];
//                const tensor coeffsD = icds[faceI];
//                const tensor coeffsL = icls[faceI];
//                const tensor coeffsU = icus[faceI];

//                diag[c](dirI + 0, dirJ + X) += coeffsD.xx() * s.x();        // xx
//                diag[c](dirI + 0, dirJ + Y) += coeffsD.xy() * s.x();        // xx
//                diag[c](dirI + 0, dirJ + Z) += coeffsD.xz() * s.x();        // xx

//                diag[c](dirI + 1, dirJ + X) += 0.5 * coeffsD.xx() * s.y();  // xy
//                diag[c](dirI + 1, dirJ + Y) += 0.5 * coeffsD.xy() * s.y();  // xy
//                diag[c](dirI + 1, dirJ + Z) += 0.5 * coeffsD.xz() * s.y();  // xy

//                diag[c](dirI + 1, dirJ + X) += 0.5 * coeffsD.yx() * s.x();  // yx
//                diag[c](dirI + 1, dirJ + Y) += 0.5 * coeffsD.yy() * s.x();  // yx
//                diag[c](dirI + 1, dirJ + Z) += 0.5 * coeffsD.yz() * s.x();  // yx

//                diag[c](dirI + 2, dirJ + X) += 0.5 * coeffsD.xx() * s.z();  // xz
//                diag[c](dirI + 2, dirJ + Y) += 0.5 * coeffsD.xy() * s.z();  // xz
//                diag[c](dirI + 2, dirJ + Z) += 0.5 * coeffsD.xz() * s.z();  // xz

//                diag[c](dirI + 2, dirJ + X) += 0.5 * coeffsD.zx() * s.x();  // zx
//                diag[c](dirI + 2, dirJ + Y) += 0.5 * coeffsD.zy() * s.x();  // zx
//                diag[c](dirI + 2, dirJ + Z) += 0.5 * coeffsD.zz() * s.x();  // zx

//                diag[c](dirI + 3, dirJ + X) += coeffsD.yx() * s.y();        // yy
//                diag[c](dirI + 3, dirJ + Y) += coeffsD.yy() * s.y();        // yy
//                diag[c](dirI + 3, dirJ + Z) += coeffsD.yz() * s.y();        // yy

//                diag[c](dirI + 4, dirJ + X) += 0.5 * coeffsD.yx() * s.z();  // yz
//                diag[c](dirI + 4, dirJ + Y) += 0.5 * coeffsD.yy() * s.z();  // yz
//                diag[c](dirI + 4, dirJ + Z) += 0.5 * coeffsD.yz() * s.z();  // yz

//                diag[c](dirI + 4, dirJ + X) += 0.5 * coeffsD.zx() * s.y();  // zy
//                diag[c](dirI + 4, dirJ + Y) += 0.5 * coeffsD.zy() * s.y();  // zy
//                diag[c](dirI + 4, dirJ + Z) += 0.5 * coeffsD.zz() * s.y();  // zy

//                diag[c](dirI + 5, dirJ + X) += coeffsD.zx() * s.z();        // zz
//                diag[c](dirI + 5, dirJ + Y) += coeffsD.zy() * s.z();        // zz
//                diag[c](dirI + 5, dirJ + Z) += coeffsD.zz() * s.z();        // zz

//                upper[secondFace](dirI + 0, dirJ + X) += coeffsU.xx() * s.x();        // xx
//                upper[secondFace](dirI + 0, dirJ + Y) += coeffsU.xy() * s.x();        // xx
//                upper[secondFace](dirI + 0, dirJ + Z) += coeffsU.xz() * s.x();        // xx

//                upper[secondFace](dirI + 1, dirJ + X) += 0.5 * coeffsU.xx() * s.y();  // xy
//                upper[secondFace](dirI + 1, dirJ + Y) += 0.5 * coeffsU.xy() * s.y();  // xy
//                upper[secondFace](dirI + 1, dirJ + Z) += 0.5 * coeffsU.xz() * s.y();  // xy

//                upper[secondFace](dirI + 1, dirJ + X) += 0.5 * coeffsU.yx() * s.x();  // yx
//                upper[secondFace](dirI + 1, dirJ + Y) += 0.5 * coeffsU.yy() * s.x();  // yx
//                upper[secondFace](dirI + 1, dirJ + Z) += 0.5 * coeffsU.yz() * s.x();  // yx

//                upper[secondFace](dirI + 2, dirJ + X) += 0.5 * coeffsU.xx() * s.z();  // xz
//                upper[secondFace](dirI + 2, dirJ + Y) += 0.5 * coeffsU.xy() * s.z();  // xz
//                upper[secondFace](dirI + 2, dirJ + Z) += 0.5 * coeffsU.xz() * s.z();  // xz

//                upper[secondFace](dirI + 2, dirJ + X) += 0.5 * coeffsU.zx() * s.x();  // zx
//                upper[secondFace](dirI + 2, dirJ + Y) += 0.5 * coeffsU.zy() * s.x();  // zx
//                upper[secondFace](dirI + 2, dirJ + Z) += 0.5 * coeffsU.zz() * s.x();  // zx

//                upper[secondFace](dirI + 3, dirJ + X) += coeffsU.yx() * s.y();        // yy
//                upper[secondFace](dirI + 3, dirJ + Y) += coeffsU.yy() * s.y();        // yy
//                upper[secondFace](dirI + 3, dirJ + Z) += coeffsU.yz() * s.y();        // yy

//                upper[secondFace](dirI + 4, dirJ + X) += 0.5 * coeffsU.yx() * s.z();  // yz
//                upper[secondFace](dirI + 4, dirJ + Y) += 0.5 * coeffsU.yy() * s.z();  // yz
//                upper[secondFace](dirI + 4, dirJ + Z) += 0.5 * coeffsU.yz() * s.z();  // yz

//                upper[secondFace](dirI + 4, dirJ + X) += 0.5 * coeffsU.zx() * s.y();  // zy
//                upper[secondFace](dirI + 4, dirJ + Y) += 0.5 * coeffsU.zy() * s.y();  // zy
//                upper[secondFace](dirI + 4, dirJ + Z) += 0.5 * coeffsU.zz() * s.y();  // zy

//                upper[secondFace](dirI + 5, dirJ + X) += coeffsU.zx() * s.z();        // zz
//                upper[secondFace](dirI + 5, dirJ + Y) += coeffsU.zy() * s.z();        // zz
//                upper[secondFace](dirI + 5, dirJ + Z) += coeffsU.zz() * s.z();        // zz

//                lower[secondFace](dirI + 0, dirJ + X) += coeffsL.xx() * s.x();        // xx
//                lower[secondFace](dirI + 0, dirJ + Y) += coeffsL.xy() * s.x();        // xx
//                lower[secondFace](dirI + 0, dirJ + Z) += coeffsL.xz() * s.x();        // xx

//                lower[secondFace](dirI + 1, dirJ + X) += 0.5 * coeffsL.xx() * s.y();  // xy
//                lower[secondFace](dirI + 1, dirJ + Y) += 0.5 * coeffsL.xy() * s.y();  // xy
//                lower[secondFace](dirI + 1, dirJ + Z) += 0.5 * coeffsL.xz() * s.y();  // xy

//                lower[secondFace](dirI + 1, dirJ + X) += 0.5 * coeffsL.yx() * s.x();  // yx
//                lower[secondFace](dirI + 1, dirJ + Y) += 0.5 * coeffsL.yy() * s.x();  // yx
//                lower[secondFace](dirI + 1, dirJ + Z) += 0.5 * coeffsL.yz() * s.x();  // yx

//                lower[secondFace](dirI + 2, dirJ + X) += 0.5 * coeffsL.xx() * s.z();  // xz
//                lower[secondFace](dirI + 2, dirJ + Y) += 0.5 * coeffsL.xy() * s.z();  // xz
//                lower[secondFace](dirI + 2, dirJ + Z) += 0.5 * coeffsL.xz() * s.z();  // xz

//                lower[secondFace](dirI + 2, dirJ + X) += 0.5 * coeffsL.zx() * s.x();  // zx
//                lower[secondFace](dirI + 2, dirJ + Y) += 0.5 * coeffsL.zy() * s.x();  // zx
//                lower[secondFace](dirI + 2, dirJ + Z) += 0.5 * coeffsL.zz() * s.x();  // zx

//                lower[secondFace](dirI + 3, dirJ + X) += coeffsL.yx() * s.y();        // yy
//                lower[secondFace](dirI + 3, dirJ + Y) += coeffsL.yy() * s.y();        // yy
//                lower[secondFace](dirI + 3, dirJ + Z) += coeffsL.yz() * s.y();        // yy

//                lower[secondFace](dirI + 4, dirJ + X) += 0.5 * coeffsL.yx() * s.z();  // yz
//                lower[secondFace](dirI + 4, dirJ + Y) += 0.5 * coeffsL.yy() * s.z();  // yz
//                lower[secondFace](dirI + 4, dirJ + Z) += 0.5 * coeffsL.yz() * s.z();  // yz

//                lower[secondFace](dirI + 4, dirJ + X) += 0.5 * coeffsL.zx() * s.y();  // zy
//                lower[secondFace](dirI + 4, dirJ + Y) += 0.5 * coeffsL.zy() * s.y();  // zy
//                lower[secondFace](dirI + 4, dirJ + Z) += 0.5 * coeffsL.zz() * s.y();  // zy

//                lower[secondFace](dirI + 5, dirJ + X) += coeffsL.zx() * s.z();        // zz
//                lower[secondFace](dirI + 5, dirJ + Y) += coeffsL.zy() * s.z();        // zz
//                lower[secondFace](dirI + 5, dirJ + Z) += coeffsL.zz() * s.z();        // zz

//            }

//        }
//        else{

//            tmp<vectorField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//            tmp<vectorField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
//            const vectorField& ic = tic();     // internal coefficient
//            const vectorField& bc = tbc();     // boundary coefficient

//            const fvPatch& patch = patchVolField.patch();   // reference to patch

//            tmp<vectorField> tsn = patch.Sf();
//            const vectorField sn = tsn();   // patch normals

//            forAll(patchVolField, faceI){   // loop all faces

//                label c = patch.faceCells()[faceI];     // boundary cell
//                const vector ci = ic[faceI];
//                const vector cb = bc[faceI];
//                const vector s = sn[faceI];
//                //xx
//                diag[c](dirI + 0, dirJ + X) += ci.x() * s.x();
//                // xy
//                diag[c](dirI + 1, dirJ + X) += 0.5 * ci.x() * s.y();
//                diag[c](dirI + 1, dirJ + Y) += 0.5 * ci.y() * s.x();
//                // xz
//                diag[c](dirI + 2, dirJ + X) += 0.5 * ci.x() * s.z();
//                diag[c](dirI + 2, dirJ + Z) += 0.5 * ci.z() * s.x();
//                // yy
//                diag[c](dirI + 3, dirJ + Y) += ci.y() * s.y();
//                // yz
//                diag[c](dirI + 4, dirJ + Y) += 0.5 * ci.y() * s.z();
//                diag[c](dirI + 4, dirJ + Z) += 0.5 * ci.z() * s.y();
//                // zz
//                diag[c](dirI + 5, dirJ + Z) += ci.z() * s.z();

//                // xx
//                B[c](dirI + 0) -= cb.x() * s.x();
//                // xy
//                B[c](dirI + 1) -= 0.5 * cb.x() * s.y();
//                B[c](dirI + 1) -= 0.5 * cb.y() * s.x();
//                // xz
//                B[c](dirI + 2) -= 0.5 * cb.x() * s.z();
//                B[c](dirI + 2) -= 0.5 * cb.z() * s.x();
//                // yy
//                B[c](dirI + 3) -= cb.y() * s.y();
//                // yz
//                B[c](dirI + 4) -= 0.5 * cb.y() * s.z();
//                B[c](dirI + 4) -= 0.5 * cb.z() * s.y();
//                // zz
//                B[c](dirI + 5) -= cb.z() * s.z();

//            }
//        }
//    }

}

template<class vectorType>
void gaussGrad(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    gaussGrad(vf,1,M,B,dirI,dirJ);


//    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

//    const word name("grad(" + vf.name() + ")");
//    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
//    const surfaceScalarField& weights = tweights();

//    const int X = 0;
//    const int Y = 1;
//    const int Z = 2;

//    // contribution of internal Field ...
//    for(int i=0; i<vf.mesh().owner().size(); i++)
//    {
//        int o = vf.mesh().owner()[i];
//        int n = vf.mesh().neighbour()[i];
//        scalar w = weights.internalField()[i];
//        vector s = vf.mesh().Sf()[i];

//        diag[o](dirI + 0, dirJ + X) += s.x() * w;
//        diag[o](dirI + 1, dirJ + X) += s.y() * w;
//        diag[o](dirI + 2, dirJ + X) += s.z() * w;
//        diag[o](dirI + 3, dirJ + Y) += s.x() * w;
//        diag[o](dirI + 4, dirJ + Y) += s.y() * w;
//        diag[o](dirI + 5, dirJ + Y) += s.z() * w;
//        diag[o](dirI + 6, dirJ + Z) += s.x() * w;
//        diag[o](dirI + 7, dirJ + Z) += s.y() * w;
//        diag[o](dirI + 8, dirJ + Z) += s.z() * w;

//        diag[n](dirI + 0, dirJ + X) -= s.x() * (1-w);
//        diag[n](dirI + 1, dirJ + X) -= s.y() * (1-w);
//        diag[n](dirI + 2, dirJ + X) -= s.z() * (1-w);
//        diag[n](dirI + 3, dirJ + Y) -= s.x() * (1-w);
//        diag[n](dirI + 4, dirJ + Y) -= s.y() * (1-w);
//        diag[n](dirI + 5, dirJ + Y) -= s.z() * (1-w);
//        diag[n](dirI + 6, dirJ + Z) -= s.x() * (1-w);
//        diag[n](dirI + 7, dirJ + Z) -= s.y() * (1-w);
//        diag[n](dirI + 8, dirJ + Z) -= s.z() * (1-w);

//        upper[i](dirI + 0, dirJ + X) += s.x() * (1-w);
//        upper[i](dirI + 1, dirJ + X) += s.y() * (1-w);
//        upper[i](dirI + 2, dirJ + X) += s.z() * (1-w);
//        upper[i](dirI + 3, dirJ + Y) += s.x() * (1-w);
//        upper[i](dirI + 4, dirJ + Y) += s.y() * (1-w);
//        upper[i](dirI + 5, dirJ + Y) += s.z() * (1-w);
//        upper[i](dirI + 6, dirJ + Z) += s.x() * (1-w);
//        upper[i](dirI + 7, dirJ + Z) += s.y() * (1-w);
//        upper[i](dirI + 8, dirJ + Z) += s.z() * (1-w);

//        lower[i](dirI + 0, dirJ + X) -= s.x() * w;
//        lower[i](dirI + 1, dirJ + X) -= s.y() * w;
//        lower[i](dirI + 2, dirJ + X) -= s.z() * w;
//        lower[i](dirI + 3, dirJ + Y) -= s.x() * w;
//        lower[i](dirI + 4, dirJ + Y) -= s.y() * w;
//        lower[i](dirI + 5, dirJ + Y) -= s.z() * w;
//        lower[i](dirI + 6, dirJ + Z) -= s.x() * w;
//        lower[i](dirI + 7, dirJ + Z) -= s.y() * w;
//        lower[i](dirI + 8, dirJ + Z) -= s.z() * w;

//    }

//    // CONTRIBUTION BOUNDARY FIELD ...
//    //
//    //      vfboundFace = ic * vfboundCell  +  bc

//    vf.boundaryField().updateCoeffs();
//    forAll(vf.boundaryField(), patchI){ // loop all patches

//        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
//        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

//        tmp<vectorField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//        tmp<vectorField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
//        const vectorField& ic = tic();     // internal coefficient
//        const vectorField& bc = tbc();     // boundary coefficient

//        const fvPatch& patch = patchVolField.patch();   // reference to patch

//        tmp<vectorField> tsn = patch.Sf();
//        const vectorField sn = tsn();   // patch normals

//        forAll(patchVolField, faceI){   // loop all faces

//            label c = patch.faceCells()[faceI];     // boundary cell
//            const vector ci = ic[faceI];
//            const vector cb = bc[faceI];
//            const vector s = sn[faceI];

//            diag[c](dirI + 0, dirJ + X) += ci.x() * s.x();
//            diag[c](dirI + 1, dirJ + X) += ci.x() * s.y();
//            diag[c](dirI + 2, dirJ + X) += ci.x() * s.z();
//            diag[c](dirI + 3, dirJ + Y) += ci.y() * s.x();
//            diag[c](dirI + 4, dirJ + Y) += ci.y() * s.y();
//            diag[c](dirI + 5, dirJ + Y) += ci.y() * s.z();
//            diag[c](dirI + 6, dirJ + Z) += ci.z() * s.x();
//            diag[c](dirI + 7, dirJ + Z) += ci.z() * s.y();
//            diag[c](dirI + 8, dirJ + Z) += ci.z() * s.z();

//            B[c](dirI + 0) -= cb.x() * s.x();
//            B[c](dirI + 1) -= cb.x() * s.y();
//            B[c](dirI + 2) -= cb.x() * s.z();
//            B[c](dirI + 3) -= cb.y() * s.x();
//            B[c](dirI + 4) -= cb.y() * s.y();
//            B[c](dirI + 5) -= cb.y() * s.z();
//            B[c](dirI + 6) -= cb.z() * s.x();
//            B[c](dirI + 7) -= cb.z() * s.y();
//            B[c](dirI + 8) -= cb.z() * s.z();

//        }

//    }
}

template<class vectorType>
void gaussGrad(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("grad(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    // contribution of internal Field ...
    for(int i=0; i<vf.mesh().owner().size(); i++)
    {
        int o = vf.mesh().owner()[i];
        int n = vf.mesh().neighbour()[i];
        scalar w = weights.internalField()[i];
        vector s = vf.mesh().Sf()[i];

        diag[o](dirI + 0, dirJ) += f*(s*w).x() * volInv[o];
        diag[o](dirI + 1, dirJ) += f*(s*w).y() * volInv[o];
        diag[o](dirI + 2, dirJ) += f*(s*w).z() * volInv[o];

        diag[n](dirI + 0, dirJ) -= f*(s*(1-w)).x() * volInv[n];
        diag[n](dirI + 1, dirJ) -= f*(s*(1-w)).y() * volInv[n];
        diag[n](dirI + 2, dirJ) -= f*(s*(1-w)).z() * volInv[n];

        upper[i](dirI + 0, dirJ) += f*(s*(1-w)).x() * volInv[o];
        upper[i](dirI + 1, dirJ) += f*(s*(1-w)).y() * volInv[o];
        upper[i](dirI + 2, dirJ) += f*(s*(1-w)).z() * volInv[o];

        lower[i](dirI + 0, dirJ) -= f*(s*w).x() * volInv[n];
        lower[i](dirI + 1, dirJ) -= f*(s*w).y() * volInv[n];
        lower[i](dirI + 2, dirJ) -= f*(s*w).z() * volInv[n];

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<scalar>& patchVolField = vf.boundaryField()[patchI];
        const fvPatch& patch = patchVolField.patch();   // reference to patch
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        if (patchVolField.type() == "implicitExtrapolation"){

            tmp<scalarField > ticd = patchVolField.valueInternalCoeffs(weightsPatchVolField);
            tmp<scalarField > ticu = patchVolField.valueInternalCoeffsUpper(weightsPatchVolField);
            tmp<scalarField > ticl = patchVolField.valueInternalCoeffsLower(weightsPatchVolField);
            const scalarField& icd = ticd();      // internal coefficient
            const scalarField& icu = ticu();
            const scalarField& icl = ticl();
            const List<label>& secondFaces = patchVolField.secondFaces();

            forAll(patchVolField, faceI){       // loop all faces

                label c = patch.faceCells()[faceI];     // boundary cell
                const label secondFace = secondFaces[faceI];

                diag[c](dirI + 0, dirJ) += f*icd[faceI]*sn[faceI].x() * volInv[c];
                diag[c](dirI + 1, dirJ) += f*icd[faceI]*sn[faceI].y() * volInv[c];
                diag[c](dirI + 2, dirJ) += f*icd[faceI]*sn[faceI].z() * volInv[c];

                upper[secondFace](dirI + 0, dirJ) += f*icu[faceI]*sn[faceI].x() * volInv[c];
                upper[secondFace](dirI + 1, dirJ) += f*icu[faceI]*sn[faceI].y() * volInv[c];
                upper[secondFace](dirI + 2, dirJ) += f*icu[faceI]*sn[faceI].z() * volInv[c];

                lower[secondFace](dirI + 0, dirJ) += f*icl[faceI]*sn[faceI].x() * volInv[c];
                lower[secondFace](dirI + 1, dirJ) += f*icl[faceI]*sn[faceI].y() * volInv[c];
                lower[secondFace](dirI + 2, dirJ) += f*icl[faceI]*sn[faceI].z() * volInv[c];

            }

        }
        else{

            tmp<scalarField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
            tmp<scalarField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
            const scalarField& ic = tic();     // internal coefficient
            const scalarField& bc = tbc();     // boundary coefficient

            forAll(patchVolField, faceI){   // loop all faces

                label c = patch.faceCells()[faceI];     // boundary cell

                diag[c](dirI + 0, dirJ) += f*ic[faceI]*sn[faceI].x() * volInv[c];
                diag[c](dirI + 1, dirJ) += f*ic[faceI]*sn[faceI].y() * volInv[c];
                diag[c](dirI + 2, dirJ) += f*ic[faceI]*sn[faceI].z() * volInv[c];

                B[c](dirI + 0) -= f*bc[faceI]*sn[faceI].x() * volInv[c];
                B[c](dirI + 1) -= f*bc[faceI]*sn[faceI].y() * volInv[c];
                B[c](dirI + 2) -= f*bc[faceI]*sn[faceI].z() * volInv[c];

            }
        }


    }

}

template<class vectorType>
void gaussSymmGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("grad(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    // contribution of internal Field ...
    for(int i=0; i<vf.mesh().owner().size(); i++)
    {
        int o = vf.mesh().owner()[i];
        int n = vf.mesh().neighbour()[i];
        scalar w = weights.internalField()[i];
        vector s = vf.mesh().Sf()[i];

        // xx
        diag[o](dirI + 0, dirJ + X) += f * s.x() * w * volInv[o];
        // xy
        diag[o](dirI + 1, dirJ + X) += f * 0.5 * s.y() * w * volInv[o];
        diag[o](dirI + 1, dirJ + Y) += f * 0.5 * s.x() * w * volInv[o];
        // xz
        diag[o](dirI + 2, dirJ + X) += f * 0.5 * s.z() * w * volInv[o];
        diag[o](dirI + 2, dirJ + Z) += f * 0.5 * s.x() * w * volInv[o];
        // yy
        diag[o](dirI + 3, dirJ + Y) += f * s.y() * w * volInv[o];
        // yz
        diag[o](dirI + 4, dirJ + Y) += f * 0.5 * s.z() * w * volInv[o];
        diag[o](dirI + 4, dirJ + Z) += f * 0.5 * s.y() * w * volInv[o];
        // zz
        diag[o](dirI + 5, dirJ + Z) += f * s.z() * w * volInv[o];

        // xx
        diag[n](dirI + 0, dirJ + X) -= f * s.x() * (1-w) * volInv[n];
        // xy
        diag[n](dirI + 1, dirJ + X) -= f * 0.5 * s.y() * (1-w) * volInv[n];
        diag[n](dirI + 1, dirJ + Y) -= f * 0.5 * s.x() * (1-w) * volInv[n];
        // xz
        diag[n](dirI + 2, dirJ + X) -= f * 0.5 * s.z() * (1-w) * volInv[n];
        diag[n](dirI + 2, dirJ + Z) -= f * 0.5 * s.x() * (1-w) * volInv[n];
        // yy
        diag[n](dirI + 3, dirJ + Y) -= f * s.y() * (1-w) * volInv[n];
        // yz
        diag[n](dirI + 4, dirJ + Y) -= f * 0.5 * s.z() * (1-w) * volInv[n];
        diag[n](dirI + 4, dirJ + Z) -= f * 0.5 * s.y() * (1-w) * volInv[n];
        // zz
        diag[n](dirI + 5, dirJ + Z) -= f * s.z() * (1-w) * volInv[n];

        // xx
        upper[i](dirI + 0, dirJ + X) += f * s.x() * (1-w) * volInv[o];
        // xy
        upper[i](dirI + 1, dirJ + X) += f * 0.5 * s.y() * (1-w) * volInv[o];
        upper[i](dirI + 1, dirJ + Y) += f * 0.5 * s.x() * (1-w) * volInv[o];
        // xz
        upper[i](dirI + 2, dirJ + X) += f * 0.5 * s.z() * (1-w) * volInv[o];
        upper[i](dirI + 2, dirJ + Z) += f * 0.5 * s.x() * (1-w) * volInv[o];
        // yy
        upper[i](dirI + 3, dirJ + Y) += f * s.y() * (1-w) * volInv[o];
        // yz
        upper[i](dirI + 4, dirJ + Y) += f * 0.5 * s.z() * (1-w) * volInv[o];
        upper[i](dirI + 4, dirJ + Z) += f * 0.5 * s.y() * (1-w) * volInv[o];
        // zz
        upper[i](dirI + 5, dirJ + Z) += f * s.z() * (1-w) * volInv[o];

        // xx
        lower[i](dirI + 0, dirJ + X) -= f * s.x() * w * volInv[n];
        // xy
        lower[i](dirI + 1, dirJ + X) -= f * 0.5 * s.y() * w * volInv[n];
        lower[i](dirI + 1, dirJ + Y) -= f * 0.5 * s.x() * w * volInv[n];
        // xz
        lower[i](dirI + 2, dirJ + X) -= f * 0.5 * s.z() * w * volInv[n];
        lower[i](dirI + 2, dirJ + Z) -= f * 0.5 * s.x() * w * volInv[n];
        // yy
        lower[i](dirI + 3, dirJ + Y) -= f * s.y() * w * volInv[n];
        // yz
        lower[i](dirI + 4, dirJ + Y) -= f * 0.5 * s.z() * w * volInv[n];
        lower[i](dirI + 4, dirJ + Z) -= f * 0.5 * s.y() * w * volInv[n];
        // zz
        lower[i](dirI + 5, dirJ + Z) -= f * s.z() * w * volInv[n];

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();

    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
        const fvPatch& patch = patchVolField.patch();   // reference to patch
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        if (patchVolField.type() == "implicitExtrapolation"){

            forAll(patchVolField, faceI){       // loop all faces

                tmp<vectorField > ticd = patchVolField.valueInternalCoeffs(weightsPatchVolField);
                tmp<vectorField > ticu = patchVolField.valueInternalCoeffsUpper(weightsPatchVolField);
                tmp<vectorField > ticl = patchVolField.valueInternalCoeffsLower(weightsPatchVolField);
                const vectorField& icd = ticd();      // internal coefficient
                const vectorField& icu = ticu();
                const vectorField& icl = ticl();
                const List<label>& secondFaces = patchVolField.secondFaces();

                label c = patch.faceCells()[faceI];     // boundary cell
                const label secondFace = secondFaces[faceI];
                const vector s = sn[faceI];

                diag[c](dirI + 0, dirJ + X) += f*icd[faceI].x() * s.x() * volInv[c];      // xx
                diag[c](dirI + 1, dirJ + X) += f*0.5 * icd[faceI].x() * s.y() * volInv[c];// xy
                diag[c](dirI + 1, dirJ + Y) += f*0.5 * icd[faceI].y() * s.x() * volInv[c];// xy
                diag[c](dirI + 2, dirJ + X) += f*0.5 * icd[faceI].x() * s.z() * volInv[c];// xz
                diag[c](dirI + 2, dirJ + Z) += f*0.5 * icd[faceI].z() * s.x() * volInv[c];// xz
                diag[c](dirI + 3, dirJ + Y) += f*icd[faceI].y() * s.y() * volInv[c];      // yy
                diag[c](dirI + 4, dirJ + Y) += f*0.5 * icd[faceI].y() * s.z() * volInv[c];// yz
                diag[c](dirI + 4, dirJ + Z) += f*0.5 * icd[faceI].z() * s.y() * volInv[c];// yz
                diag[c](dirI + 5, dirJ + Z) += f*icd[faceI].z() * s.z() * volInv[c];      // zz

                if(icu[faceI] != vector(0,0,0)){
                    upper[secondFace](dirI + 0, dirJ + X) += f*icu[faceI].x() * s.x() * volInv[c];      // xx
                    upper[secondFace](dirI + 1, dirJ + X) += f*0.5 * icu[faceI].x() * s.y() * volInv[c];// xy
                    upper[secondFace](dirI + 1, dirJ + Y) += f*0.5 * icu[faceI].y() * s.x() * volInv[c];// xy
                    upper[secondFace](dirI + 2, dirJ + X) += f*0.5 * icu[faceI].x() * s.z() * volInv[c];// xz
                    upper[secondFace](dirI + 2, dirJ + Z) += f*0.5 * icu[faceI].z() * s.x() * volInv[c];// xz
                    upper[secondFace](dirI + 3, dirJ + Y) += f*icu[faceI].y() * s.y() * volInv[c];      // yy
                    upper[secondFace](dirI + 4, dirJ + Y) += f*0.5 * icu[faceI].y() * s.z() * volInv[c];// yz
                    upper[secondFace](dirI + 4, dirJ + Z) += f*0.5 * icu[faceI].z() * s.y() * volInv[c];// yz
                    upper[secondFace](dirI + 5, dirJ + Z) += f*icu[faceI].z() * s.z() * volInv[c];      // zz
                }
                if(icl[faceI] != vector(0,0,0)){
                    lower[secondFace](dirI + 0, dirJ + X) += f*icl[faceI].x() * s.x() * volInv[c];      // xx
                    lower[secondFace](dirI + 1, dirJ + X) += f*0.5 * icl[faceI].x() * s.y() * volInv[c];// xy
                    lower[secondFace](dirI + 1, dirJ + Y) += f*0.5 * icl[faceI].y() * s.x() * volInv[c];// xy
                    lower[secondFace](dirI + 2, dirJ + X) += f*0.5 * icl[faceI].x() * s.z() * volInv[c];// xz
                    lower[secondFace](dirI + 2, dirJ + Z) += f*0.5 * icl[faceI].z() * s.x() * volInv[c];// xz
                    lower[secondFace](dirI + 3, dirJ + Y) += f*icl[faceI].y() * s.y() * volInv[c];      // yy
                    lower[secondFace](dirI + 4, dirJ + Y) += f*0.5 * icl[faceI].y() * s.z() * volInv[c];// yz
                    lower[secondFace](dirI + 4, dirJ + Z) += f*0.5 * icl[faceI].z() * s.y() * volInv[c];// yz
                    lower[secondFace](dirI + 5, dirJ + Z) += f*icl[faceI].z() * s.z() * volInv[c];      // zz
                }

            }

        }
        else if (patchVolField.type() == "myImplicitVelocity"){

            tmp<tensorField > ticds = patchVolField.valueInternalCoeffsTensor(weightsPatchVolField);
            tmp<tensorField > ticus = patchVolField.valueInternalCoeffsUpperTensor(weightsPatchVolField);
            tmp<tensorField > ticls = patchVolField.valueInternalCoeffsLowerTensor(weightsPatchVolField);
            const tensorField& icds = ticds();      // internal coefficient
            const tensorField& icus = ticus();
            const tensorField& icls = ticls();
            const List<label>& secondFaces = patchVolField.secondFaces();

            forAll(patchVolField, faceI){       // loop all faces

                label c = patch.faceCells()[faceI];     // boundary cell
                const label secondFace = secondFaces[faceI];
                const vector s = sn[faceI];
                const tensor coeffsD = icds[faceI];
                const tensor coeffsL = icls[faceI];
                const tensor coeffsU = icus[faceI];

                diag[c](dirI + 0, dirJ + X) += f*coeffsD.xx() * s.x() * volInv[c];        // xx
                diag[c](dirI + 0, dirJ + Y) += f*coeffsD.xy() * s.x() * volInv[c];        // xx
                diag[c](dirI + 0, dirJ + Z) += f*coeffsD.xz() * s.x() * volInv[c];        // xx

                diag[c](dirI + 1, dirJ + X) += f*0.5 * coeffsD.xx() * s.y() * volInv[c];  // xy
                diag[c](dirI + 1, dirJ + Y) += f*0.5 * coeffsD.xy() * s.y() * volInv[c];  // xy
                diag[c](dirI + 1, dirJ + Z) += f*0.5 * coeffsD.xz() * s.y() * volInv[c];  // xy

                diag[c](dirI + 1, dirJ + X) += f*0.5 * coeffsD.yx() * s.x() * volInv[c];  // yx
                diag[c](dirI + 1, dirJ + Y) += f*0.5 * coeffsD.yy() * s.x() * volInv[c];  // yx
                diag[c](dirI + 1, dirJ + Z) += f*0.5 * coeffsD.yz() * s.x() * volInv[c];  // yx

                diag[c](dirI + 2, dirJ + X) += f*0.5 * coeffsD.xx() * s.z() * volInv[c];  // xz
                diag[c](dirI + 2, dirJ + Y) += f*0.5 * coeffsD.xy() * s.z() * volInv[c];  // xz
                diag[c](dirI + 2, dirJ + Z) += f*0.5 * coeffsD.xz() * s.z() * volInv[c];  // xz

                diag[c](dirI + 2, dirJ + X) += f*0.5 * coeffsD.zx() * s.x() * volInv[c];  // zx
                diag[c](dirI + 2, dirJ + Y) += f*0.5 * coeffsD.zy() * s.x() * volInv[c];  // zx
                diag[c](dirI + 2, dirJ + Z) += f*0.5 * coeffsD.zz() * s.x() * volInv[c];  // zx

                diag[c](dirI + 3, dirJ + X) += f*coeffsD.yx() * s.y() * volInv[c];        // yy
                diag[c](dirI + 3, dirJ + Y) += f*coeffsD.yy() * s.y() * volInv[c];        // yy
                diag[c](dirI + 3, dirJ + Z) += f*coeffsD.yz() * s.y() * volInv[c];        // yy

                diag[c](dirI + 4, dirJ + X) += f*0.5 * coeffsD.yx() * s.z() * volInv[c];  // yz
                diag[c](dirI + 4, dirJ + Y) += f*0.5 * coeffsD.yy() * s.z() * volInv[c];  // yz
                diag[c](dirI + 4, dirJ + Z) += f*0.5 * coeffsD.yz() * s.z() * volInv[c];  // yz

                diag[c](dirI + 4, dirJ + X) += f*0.5 * coeffsD.zx() * s.y() * volInv[c];  // zy
                diag[c](dirI + 4, dirJ + Y) += f*0.5 * coeffsD.zy() * s.y() * volInv[c];  // zy
                diag[c](dirI + 4, dirJ + Z) += f*0.5 * coeffsD.zz() * s.y() * volInv[c];  // zy

                diag[c](dirI + 5, dirJ + X) += f*coeffsD.zx() * s.z() * volInv[c];        // zz
                diag[c](dirI + 5, dirJ + Y) += f*coeffsD.zy() * s.z() * volInv[c];        // zz
                diag[c](dirI + 5, dirJ + Z) += f*coeffsD.zz() * s.z() * volInv[c];        // zz

                upper[secondFace](dirI + 0, dirJ + X) += f*coeffsU.xx() * s.x() * volInv[c];        // xx
                upper[secondFace](dirI + 0, dirJ + Y) += f*coeffsU.xy() * s.x() * volInv[c];        // xx
                upper[secondFace](dirI + 0, dirJ + Z) += f*coeffsU.xz() * s.x() * volInv[c];        // xx

                upper[secondFace](dirI + 1, dirJ + X) += f*0.5 * coeffsU.xx() * s.y() * volInv[c];  // xy
                upper[secondFace](dirI + 1, dirJ + Y) += f*0.5 * coeffsU.xy() * s.y() * volInv[c];  // xy
                upper[secondFace](dirI + 1, dirJ + Z) += f*0.5 * coeffsU.xz() * s.y() * volInv[c];  // xy

                upper[secondFace](dirI + 1, dirJ + X) += f*0.5 * coeffsU.yx() * s.x() * volInv[c];  // yx
                upper[secondFace](dirI + 1, dirJ + Y) += f*0.5 * coeffsU.yy() * s.x() * volInv[c];  // yx
                upper[secondFace](dirI + 1, dirJ + Z) += f*0.5 * coeffsU.yz() * s.x() * volInv[c];  // yx

                upper[secondFace](dirI + 2, dirJ + X) += f*0.5 * coeffsU.xx() * s.z() * volInv[c];  // xz
                upper[secondFace](dirI + 2, dirJ + Y) += f*0.5 * coeffsU.xy() * s.z() * volInv[c];  // xz
                upper[secondFace](dirI + 2, dirJ + Z) += f*0.5 * coeffsU.xz() * s.z() * volInv[c];  // xz

                upper[secondFace](dirI + 2, dirJ + X) += f*0.5 * coeffsU.zx() * s.x() * volInv[c];  // zx
                upper[secondFace](dirI + 2, dirJ + Y) += f*0.5 * coeffsU.zy() * s.x() * volInv[c];  // zx
                upper[secondFace](dirI + 2, dirJ + Z) += f*0.5 * coeffsU.zz() * s.x() * volInv[c];  // zx

                upper[secondFace](dirI + 3, dirJ + X) += f*coeffsU.yx() * s.y() * volInv[c];        // yy
                upper[secondFace](dirI + 3, dirJ + Y) += f*coeffsU.yy() * s.y() * volInv[c];        // yy
                upper[secondFace](dirI + 3, dirJ + Z) += f*coeffsU.yz() * s.y() * volInv[c];        // yy

                upper[secondFace](dirI + 4, dirJ + X) += f*0.5 * coeffsU.yx() * s.z() * volInv[c];  // yz
                upper[secondFace](dirI + 4, dirJ + Y) += f*0.5 * coeffsU.yy() * s.z() * volInv[c];  // yz
                upper[secondFace](dirI + 4, dirJ + Z) += f*0.5 * coeffsU.yz() * s.z() * volInv[c];  // yz

                upper[secondFace](dirI + 4, dirJ + X) += f*0.5 * coeffsU.zx() * s.y() * volInv[c];  // zy
                upper[secondFace](dirI + 4, dirJ + Y) += f*0.5 * coeffsU.zy() * s.y() * volInv[c];  // zy
                upper[secondFace](dirI + 4, dirJ + Z) += f*0.5 * coeffsU.zz() * s.y() * volInv[c];  // zy

                upper[secondFace](dirI + 5, dirJ + X) += f*coeffsU.zx() * s.z() * volInv[c];        // zz
                upper[secondFace](dirI + 5, dirJ + Y) += f*coeffsU.zy() * s.z() * volInv[c];        // zz
                upper[secondFace](dirI + 5, dirJ + Z) += f*coeffsU.zz() * s.z() * volInv[c];        // zz

                lower[secondFace](dirI + 0, dirJ + X) += f*coeffsL.xx() * s.x() * volInv[c];        // xx
                lower[secondFace](dirI + 0, dirJ + Y) += f*coeffsL.xy() * s.x() * volInv[c];        // xx
                lower[secondFace](dirI + 0, dirJ + Z) += f*coeffsL.xz() * s.x() * volInv[c];        // xx

                lower[secondFace](dirI + 1, dirJ + X) += f*0.5 * coeffsL.xx() * s.y() * volInv[c];  // xy
                lower[secondFace](dirI + 1, dirJ + Y) += f*0.5 * coeffsL.xy() * s.y() * volInv[c];  // xy
                lower[secondFace](dirI + 1, dirJ + Z) += f*0.5 * coeffsL.xz() * s.y() * volInv[c];  // xy

                lower[secondFace](dirI + 1, dirJ + X) += f*0.5 * coeffsL.yx() * s.x() * volInv[c];  // yx
                lower[secondFace](dirI + 1, dirJ + Y) += f*0.5 * coeffsL.yy() * s.x() * volInv[c];  // yx
                lower[secondFace](dirI + 1, dirJ + Z) += f*0.5 * coeffsL.yz() * s.x() * volInv[c];  // yx

                lower[secondFace](dirI + 2, dirJ + X) += f*0.5 * coeffsL.xx() * s.z() * volInv[c];  // xz
                lower[secondFace](dirI + 2, dirJ + Y) += f*0.5 * coeffsL.xy() * s.z() * volInv[c];  // xz
                lower[secondFace](dirI + 2, dirJ + Z) += f*0.5 * coeffsL.xz() * s.z() * volInv[c];  // xz

                lower[secondFace](dirI + 2, dirJ + X) += f*0.5 * coeffsL.zx() * s.x() * volInv[c];  // zx
                lower[secondFace](dirI + 2, dirJ + Y) += f*0.5 * coeffsL.zy() * s.x() * volInv[c];  // zx
                lower[secondFace](dirI + 2, dirJ + Z) += f*0.5 * coeffsL.zz() * s.x() * volInv[c];  // zx

                lower[secondFace](dirI + 3, dirJ + X) += f*coeffsL.yx() * s.y() * volInv[c];        // yy
                lower[secondFace](dirI + 3, dirJ + Y) += f*coeffsL.yy() * s.y() * volInv[c];        // yy
                lower[secondFace](dirI + 3, dirJ + Z) += f*coeffsL.yz() * s.y() * volInv[c];        // yy

                lower[secondFace](dirI + 4, dirJ + X) += f*0.5 * coeffsL.yx() * s.z() * volInv[c];  // yz
                lower[secondFace](dirI + 4, dirJ + Y) += f*0.5 * coeffsL.yy() * s.z() * volInv[c];  // yz
                lower[secondFace](dirI + 4, dirJ + Z) += f*0.5 * coeffsL.yz() * s.z() * volInv[c];  // yz

                lower[secondFace](dirI + 4, dirJ + X) += f*0.5 * coeffsL.zx() * s.y() * volInv[c];  // zy
                lower[secondFace](dirI + 4, dirJ + Y) += f*0.5 * coeffsL.zy() * s.y() * volInv[c];  // zy
                lower[secondFace](dirI + 4, dirJ + Z) += f*0.5 * coeffsL.zz() * s.y() * volInv[c];  // zy

                lower[secondFace](dirI + 5, dirJ + X) += f*coeffsL.zx() * s.z() * volInv[c];        // zz
                lower[secondFace](dirI + 5, dirJ + Y) += f*coeffsL.zy() * s.z() * volInv[c];        // zz
                lower[secondFace](dirI + 5, dirJ + Z) += f*coeffsL.zz() * s.z() * volInv[c];        // zz

            }

        }
        else{

            tmp<vectorField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
            tmp<vectorField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
            const vectorField& ic = tic();     // internal coefficient
            const vectorField& bc = tbc();     // boundary coefficient

            const fvPatch& patch = patchVolField.patch();   // reference to patch

            tmp<vectorField> tsn = patch.Sf();
            const vectorField sn = tsn();   // patch normals

            forAll(patchVolField, faceI){   // loop all faces

                label c = patch.faceCells()[faceI];     // boundary cell
                const vector ci = ic[faceI];
                const vector cb = bc[faceI];
                const vector s = sn[faceI];
                //xx
                diag[c](dirI + 0, dirJ + X) += f * ci.x() * s.x() * volInv[c];
                // xy
                diag[c](dirI + 1, dirJ + X) += f * 0.5 * ci.x() * s.y() * volInv[c];
                diag[c](dirI + 1, dirJ + Y) += f * 0.5 * ci.y() * s.x() * volInv[c];
                // xz
                diag[c](dirI + 2, dirJ + X) += f * 0.5 * ci.x() * s.z() * volInv[c];
                diag[c](dirI + 2, dirJ + Z) += f * 0.5 * ci.z() * s.x() * volInv[c];
                // yy
                diag[c](dirI + 3, dirJ + Y) += f * ci.y() * s.y() * volInv[c];
                // yz
                diag[c](dirI + 4, dirJ + Y) += f * 0.5 * ci.y() * s.z() * volInv[c];
                diag[c](dirI + 4, dirJ + Z) += f * 0.5 * ci.z() * s.y() * volInv[c];
                // zz
                diag[c](dirI + 5, dirJ + Z) += f * ci.z() * s.z() * volInv[c];

                // xx
                B[c](dirI + 0) -= f * cb.x() * s.x() * volInv[c];
                // xy
                B[c](dirI + 1) -= f * 0.5 * cb.x() * s.y() * volInv[c];
                B[c](dirI + 1) -= f * 0.5 * cb.y() * s.x() * volInv[c];
                // xz
                B[c](dirI + 2) -= f * 0.5 * cb.x() * s.z() * volInv[c];
                B[c](dirI + 2) -= f * 0.5 * cb.z() * s.x() * volInv[c];
                // yy
                B[c](dirI + 3) -= f * cb.y() * s.y() * volInv[c];
                // yz
                B[c](dirI + 4) -= f * 0.5 * cb.y() * s.z() * volInv[c];
                B[c](dirI + 4) -= f * 0.5 * cb.z() * s.y() * volInv[c];
                // zz
                B[c](dirI + 5) -= f * cb.z() * s.z() * volInv[c];

            }
        }
    }

}

template<class vectorType>
void gaussGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("grad(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    const int X = 0;
    const int Y = 1;
    const int Z = 2;

    // contribution of internal Field ...
    for(int i=0; i<vf.mesh().owner().size(); i++)
    {
        int o = vf.mesh().owner()[i];
        int n = vf.mesh().neighbour()[i];
        scalar w = weights.internalField()[i];
        vector s = vf.mesh().Sf()[i];

        diag[o](dirI + 0, dirJ + X) += f * s.x() * w * volInv[o];
        diag[o](dirI + 1, dirJ + X) += f * s.y() * w * volInv[o];
        diag[o](dirI + 2, dirJ + X) += f * s.z() * w * volInv[o];
        diag[o](dirI + 3, dirJ + Y) += f * s.x() * w * volInv[o];
        diag[o](dirI + 4, dirJ + Y) += f * s.y() * w * volInv[o];
        diag[o](dirI + 5, dirJ + Y) += f * s.z() * w * volInv[o];
        diag[o](dirI + 6, dirJ + Z) += f * s.x() * w * volInv[o];
        diag[o](dirI + 7, dirJ + Z) += f * s.y() * w * volInv[o];
        diag[o](dirI + 8, dirJ + Z) += f * s.z() * w * volInv[o];

        diag[n](dirI + 0, dirJ + X) -= f * s.x() * (1-w) * volInv[n];
        diag[n](dirI + 1, dirJ + X) -= f * s.y() * (1-w) * volInv[n];
        diag[n](dirI + 2, dirJ + X) -= f * s.z() * (1-w) * volInv[n];
        diag[n](dirI + 3, dirJ + Y) -= f * s.x() * (1-w) * volInv[n];
        diag[n](dirI + 4, dirJ + Y) -= f * s.y() * (1-w) * volInv[n];
        diag[n](dirI + 5, dirJ + Y) -= f * s.z() * (1-w) * volInv[n];
        diag[n](dirI + 6, dirJ + Z) -= f * s.x() * (1-w) * volInv[n];
        diag[n](dirI + 7, dirJ + Z) -= f * s.y() * (1-w) * volInv[n];
        diag[n](dirI + 8, dirJ + Z) -= f * s.z() * (1-w) * volInv[n];

        upper[i](dirI + 0, dirJ + X) += f * s.x() * (1-w) * volInv[o];
        upper[i](dirI + 1, dirJ + X) += f * s.y() * (1-w) * volInv[o];
        upper[i](dirI + 2, dirJ + X) += f * s.z() * (1-w) * volInv[o];
        upper[i](dirI + 3, dirJ + Y) += f * s.x() * (1-w) * volInv[o];
        upper[i](dirI + 4, dirJ + Y) += f * s.y() * (1-w) * volInv[o];
        upper[i](dirI + 5, dirJ + Y) += f * s.z() * (1-w) * volInv[o];
        upper[i](dirI + 6, dirJ + Z) += f * s.x() * (1-w) * volInv[o];
        upper[i](dirI + 7, dirJ + Z) += f * s.y() * (1-w) * volInv[o];
        upper[i](dirI + 8, dirJ + Z) += f * s.z() * (1-w) * volInv[o];

        lower[i](dirI + 0, dirJ + X) -= f * s.x() * w * volInv[n];
        lower[i](dirI + 1, dirJ + X) -= f * s.y() * w * volInv[n];
        lower[i](dirI + 2, dirJ + X) -= f * s.z() * w * volInv[n];
        lower[i](dirI + 3, dirJ + Y) -= f * s.x() * w * volInv[n];
        lower[i](dirI + 4, dirJ + Y) -= f * s.y() * w * volInv[n];
        lower[i](dirI + 5, dirJ + Y) -= f * s.z() * w * volInv[n];
        lower[i](dirI + 6, dirJ + Z) -= f * s.x() * w * volInv[n];
        lower[i](dirI + 7, dirJ + Z) -= f * s.y() * w * volInv[n];
        lower[i](dirI + 8, dirJ + Z) -= f * s.z() * w * volInv[n];

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with (uniform 1)

        tmp<vectorField > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<vectorField > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const vectorField& ic = tic();     // internal coefficient
        const vectorField& bc = tbc();     // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        tmp<vectorField> tsn = patch.Sf();
        const vectorField sn = tsn();   // patch normals

        forAll(patchVolField, faceI){   // loop all faces

            label c = patch.faceCells()[faceI];     // boundary cell
            const vector ci = ic[faceI];
            const vector cb = bc[faceI];
            const vector s = sn[faceI];

            diag[c](dirI + 0, dirJ + X) += f * ci.x() * s.x() * volInv[c];
            diag[c](dirI + 1, dirJ + X) += f * ci.x() * s.y() * volInv[c];
            diag[c](dirI + 2, dirJ + X) += f * ci.x() * s.z() * volInv[c];
            diag[c](dirI + 3, dirJ + Y) += f * ci.y() * s.x() * volInv[c];
            diag[c](dirI + 4, dirJ + Y) += f * ci.y() * s.y() * volInv[c];
            diag[c](dirI + 5, dirJ + Y) += f * ci.y() * s.z() * volInv[c];
            diag[c](dirI + 6, dirJ + Z) += f * ci.z() * s.x() * volInv[c];
            diag[c](dirI + 7, dirJ + Z) += f * ci.z() * s.y() * volInv[c];
            diag[c](dirI + 8, dirJ + Z) += f * ci.z() * s.z() * volInv[c];

            B[c](dirI + 0) -= f * cb.x() * s.x() * volInv[c];
            B[c](dirI + 1) -= f * cb.x() * s.y() * volInv[c];
            B[c](dirI + 2) -= f * cb.x() * s.z() * volInv[c];
            B[c](dirI + 3) -= f * cb.y() * s.x() * volInv[c];
            B[c](dirI + 4) -= f * cb.y() * s.y() * volInv[c];
            B[c](dirI + 5) -= f * cb.y() * s.z() * volInv[c];
            B[c](dirI + 6) -= f * cb.z() * s.x() * volInv[c];
            B[c](dirI + 7) -= f * cb.z() * s.y() * volInv[c];
            B[c](dirI + 8) -= f * cb.z() * s.z() * volInv[c];

        }

    }
}

}

