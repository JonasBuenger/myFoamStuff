#include"fvCFD.H"
#include "../useful/interpolationWeights.H"

namespace myFvm{

template<class vectorType>
void gaussGrad(volScalarField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("grad(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    if (!updateOnlyRHS){

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

        if (patchVolField.type() == "implicitExtrapolation" && !updateOnlyRHS){

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

                if (!updateOnlyRHS){

                    diag[c](dirI + 0, dirJ) += f*ic[faceI]*sn[faceI].x() * volInv[c];
                    diag[c](dirI + 1, dirJ) += f*ic[faceI]*sn[faceI].y() * volInv[c];
                    diag[c](dirI + 2, dirJ) += f*ic[faceI]*sn[faceI].z() * volInv[c];

                }

                B[c](dirI + 0) -= f*bc[faceI]*sn[faceI].x() * volInv[c];
                B[c](dirI + 1) -= f*bc[faceI]*sn[faceI].y() * volInv[c];
                B[c](dirI + 2) -= f*bc[faceI]*sn[faceI].z() * volInv[c];

            }
        }


    }

}

template<class vectorType>
void gaussSymmGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

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

    if (!updateOnlyRHS){

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

        if (patchVolField.type() == "implicitExtrapolation" && !updateOnlyRHS){

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
        else if (patchVolField.type() == "myImplicitVelocity" && !updateOnlyRHS){

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

                if(!updateOnlyRHS){

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

                }

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
void gaussGrad(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

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
    if (!updateOnlyRHS){

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

            if (!updateOnlyRHS){

                diag[c](dirI + 0, dirJ + X) += f * ci.x() * s.x() * volInv[c];
                diag[c](dirI + 1, dirJ + X) += f * ci.x() * s.y() * volInv[c];
                diag[c](dirI + 2, dirJ + X) += f * ci.x() * s.z() * volInv[c];
                diag[c](dirI + 3, dirJ + Y) += f * ci.y() * s.x() * volInv[c];
                diag[c](dirI + 4, dirJ + Y) += f * ci.y() * s.y() * volInv[c];
                diag[c](dirI + 5, dirJ + Y) += f * ci.y() * s.z() * volInv[c];
                diag[c](dirI + 6, dirJ + Z) += f * ci.z() * s.x() * volInv[c];
                diag[c](dirI + 7, dirJ + Z) += f * ci.z() * s.y() * volInv[c];
                diag[c](dirI + 8, dirJ + Z) += f * ci.z() * s.z() * volInv[c];

            }

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

