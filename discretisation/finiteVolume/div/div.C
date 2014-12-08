#include"fvCFD.H"
#include "../useful/interpolationWeights.H"

namespace myFvm{

template<class vectorType>
void gaussDiv(volVectorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    gaussDiv(vf,1,M,B,dirI,dirJ);


//    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

//    const word name("div(" + vf.name() + ")");
//    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
//    const surfaceScalarField& weights = tweights();
//    const volScalarField volInv = 1/vf.mesh().V();

//    for(int i=0; i<vf.mesh().owner().size(); i++){

//        int o = vf.mesh().owner()[i];
//        int n = vf.mesh().neighbour()[i];
//        scalar w = weights.internalField()[i];
//        vector s = vf.mesh().Sf()[i];

//        diag[o](dirI + 0, dirJ + 0) += w*s.x() * volInv[o];
//        diag[o](dirI + 0, dirJ + 1) += w*s.y() * volInv[o];
//        diag[o](dirI + 0, dirJ + 2) += w*s.z() * volInv[o];

//        diag[n](dirI + 0, dirJ + 0) -= (1-w)*s.x() * volInv[n];
//        diag[n](dirI + 0, dirJ + 1) -= (1-w)*s.y() * volInv[n];
//        diag[n](dirI + 0, dirJ + 2) -= (1-w)*s.z() * volInv[n];

//        upper[i](dirI + 0, dirJ + 0) += (1-w)*s.x() * volInv[o];
//        upper[i](dirI + 0, dirJ + 1) += (1-w)*s.y() * volInv[o];
//        upper[i](dirI + 0, dirJ + 2) += (1-w)*s.z() * volInv[o];

//        lower[i](dirI + 0, dirJ + 0) -= w*s.x() * volInv[n];
//        lower[i](dirI + 0, dirJ + 1) -= w*s.y() * volInv[n];
//        lower[i](dirI + 0, dirJ + 2) -= w*s.z() * volInv[n];

//    }

//    // CONTRIBUTION BOUNDARY FIELD ...
//    //
//    //      vfboundFace = ic * vfboundCell  +  bc

//    vf.boundaryField().updateCoeffs();
//    forAll(vf.boundaryField(), patchI){ // loop all patches

//        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
//        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

//        if (patchVolField.type() == "myImplicitVelocity"){
//            notImplemented("myImplicitVelocity für div(vectorField) noch nicht implementiert!!!");

//            /*
//            tmp<tensorField> ticd = patchVolField.valueInternalCoeffsTensor(weightsPatchVolField);
//            tmp<tensorField> ticu = patchVolField.valueInternalCoeffsUpperTensor(weightsPatchVolField);
//            tmp<tensorField> ticl = patchVolField.valueInternalCoeffsLowerTensor(weightsPatchVolField);
//            const tensorField& icd = ticd();              // internal coefficient
//            const tensorField& icu = ticu();              // upper coefficient
//            const tensorField& icl = ticl();              // lower coefficient

//            const labelList& faceCells = patchVolField.patch().faceCells();
//            const vectorField& sn = patchVolField.patch().nf();
//            const labelList& secondFaces = patchVolField.secondFaces();

//            forAll(patchVolField, faceI){               // loop all faces

//                const label c = faceCells[faceI];
//                const label secondFace = secondFaces[faceI];
//                const vector n = sn[faceI];
//                const tensor coeffsD = icd[faceI];
//                const tensor coeffsU = icu[faceI];
//                const tensor coeffsL = icl[faceI];


//                diag[c](dirI + 0, dirJ + 0) += (coeffsD.xx()*n.x() + coeffsD.yx()*n.y() + coeffsD.zx()*n.z()) * volInv[c];
//                diag[c](dirI + 0, dirJ + 1) += (coeffsD.xy()*n.x() + coeffsD.yy()*n.y() + coeffsD.zy()*n.z()) * volInv[c];
//                diag[c](dirI + 0, dirJ + 2) += (coeffsD.xz()*n.x() + coeffsD.yz()*n.y() + coeffsD.zz()*n.z()) * volInv[c];

//                upper[secondFace](dirI + 0, dirJ + 0) += (coeffsU.xx()*n.x() + coeffsU.yx()*n.y() + coeffsU.zx()*n.z()) * volInv[c];
//                upper[secondFace](dirI + 0, dirJ + 1) += (coeffsU.xy()*n.x() + coeffsU.yy()*n.y() + coeffsU.zy()*n.z()) * volInv[c];
//                upper[secondFace](dirI + 0, dirJ + 2) += (coeffsU.xz()*n.x() + coeffsU.yz()*n.y() + coeffsU.zz()*n.z()) * volInv[c];

//                lower[secondFace](dirI + 0, dirJ + 0) += (coeffsL.xx()*n.x() + coeffsL.yx()*n.y() + coeffsL.zx()*n.z()) * volInv[c];
//                lower[secondFace](dirI + 0, dirJ + 1) += (coeffsL.xy()*n.x() + coeffsL.yy()*n.y() + coeffsL.zy()*n.z()) * volInv[c];
//                lower[secondFace](dirI + 0, dirJ + 2) += (coeffsL.xz()*n.x() + coeffsL.yz()*n.y() + coeffsL.zz()*n.z()) * volInv[c];


//            }

//            //Hier noch füllen!!!
//        */
//        }
//        else if(patchVolField.type() == "implicitExtrapolation"){
//            notImplemented("implicitExtrapolation für div(vectorField) noch nicht implementiert!!!");
//        }
//        else{

//            tmp<vectorField> tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//            tmp<vectorField> tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
//            const vectorField& ic = tic();              // internal coefficient
//            const vectorField& bc = tbc();              // boundary coefficient

//            const fvPatch& patch = patchVolField.patch();   // reference to patch

//            const vectorField& sn = patch.Sf();         // patch normals

//            forAll(patchVolField, faceI){               // loop all faces

//                const label c = patch.faceCells()[faceI];
//                const vector v = cmptMultiply(ic[faceI],sn[faceI]);

//                diag[c](dirI + 0, dirJ + 0) += v.x() * volInv[c];
//                diag[c](dirI + 0, dirJ + 1) += v.y() * volInv[c];
//                diag[c](dirI + 0, dirJ + 2) += v.z() * volInv[c];

//                B[c](dirI) -= (bc[faceI] & sn[faceI]) * volInv[c];
//            }

//        }
//    }

}

template<class vectorType>
void gaussDiv(volSymmTensorField& vf, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    gaussDiv(vf,1,M,B,dirI,dirJ);

//    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
//    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

//    const word name("div(" + vf.name() + ")");
//    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
//    const surfaceScalarField& weights = tweights();

//    const int XX = 0;
//    const int XY = 1;
//    const int XZ = 2;
//    const int YY = 3;
//    const int YZ = 4;
//    const int ZZ = 5;

//    const int YX = XY;
//    const int ZX = XZ;
//    const int ZY = YZ;

//    for(int i=0; i<vf.mesh().owner().size(); i++){

//        int o = vf.mesh().owner()[i];
//        int n = vf.mesh().neighbour()[i];
//        scalar w = weights.internalField()[i];
//        vector s = vf.mesh().Sf()[i];

//        // zeile 1
//        diag[o](dirI + 0, dirJ + XX) += w*s.x();
//        diag[o](dirI + 0, dirJ + XY) += w*s.y();
//        diag[o](dirI + 0, dirJ + XZ) += w*s.z();

//        diag[n](dirI + 0, dirJ + XX) -= (1-w)*s.x();
//        diag[n](dirI + 0, dirJ + XY) -= (1-w)*s.y();
//        diag[n](dirI + 0, dirJ + XZ) -= (1-w)*s.z();

//        upper[i](dirI + 0, dirJ + XX) += (1-w)*s.x();
//        upper[i](dirI + 0, dirJ + XY) += (1-w)*s.y();
//        upper[i](dirI + 0, dirJ + XZ) += (1-w)*s.z();

//        lower[i](dirI + 0, dirJ + XX) -= w*s.x();
//        lower[i](dirI + 0, dirJ + XY) -= w*s.y();
//        lower[i](dirI + 0, dirJ + XZ) -= w*s.z();

//        // zeile 2
//        diag[o](dirI + 1, dirJ + YX) += w*s.x();
//        diag[o](dirI + 1, dirJ + YY) += w*s.y();
//        diag[o](dirI + 1, dirJ + YZ) += w*s.z();

//        diag[n](dirI + 1, dirJ + YX) -= (1-w)*s.x();
//        diag[n](dirI + 1, dirJ + YY) -= (1-w)*s.y();
//        diag[n](dirI + 1, dirJ + YZ) -= (1-w)*s.z();

//        upper[i](dirI + 1, dirJ + YX) += (1-w)*s.x();
//        upper[i](dirI + 1, dirJ + YY) += (1-w)*s.y();
//        upper[i](dirI + 1, dirJ + YZ) += (1-w)*s.z();

//        lower[i](dirI + 1, dirJ + YX) -= w*s.x();
//        lower[i](dirI + 1, dirJ + YY) -= w*s.y();
//        lower[i](dirI + 1, dirJ + YZ) -= w*s.z();

//        // zeile 3
//        diag[o](dirI + 2, dirJ + ZX) += w*s.x();
//        diag[o](dirI + 2, dirJ + ZY) += w*s.y();
//        diag[o](dirI + 2, dirJ + ZZ) += w*s.z();

//        diag[n](dirI + 2, dirJ + ZX) -= (1-w)*s.x();
//        diag[n](dirI + 2, dirJ + ZY) -= (1-w)*s.y();
//        diag[n](dirI + 2, dirJ + ZZ) -= (1-w)*s.z();

//        upper[i](dirI + 2, dirJ + ZX) += (1-w)*s.x();
//        upper[i](dirI + 2, dirJ + ZY) += (1-w)*s.y();
//        upper[i](dirI + 2, dirJ + ZZ) += (1-w)*s.z();

//        lower[i](dirI + 2, dirJ + ZX) -= w*s.x();
//        lower[i](dirI + 2, dirJ + ZY) -= w*s.y();
//        lower[i](dirI + 2, dirJ + ZZ) -= w*s.z();
//    }

//    // CONTRIBUTION BOUNDARY FIELD ...
//    //
//    //      vfboundFace = ic * vfboundCell  +  bc

//    vf.boundaryField().updateCoeffs();
//    forAll(vf.boundaryField(), patchI){ // loop all patches

//        const fvPatchField<symmTensor>& patchVolField = vf.boundaryField()[patchI];
//        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

//        tmp<symmTensorField> tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
//        tmp<symmTensorField> tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
//        const symmTensorField& ic = tic();              // internal coefficient
//        const symmTensorField& bc = tbc();              // boundary coefficient

//        const fvPatch& patch = patchVolField.patch();   // reference to patch

//        const vectorField& sn = patch.Sf();         // patch normals

//        forAll(patchVolField, faceI){               // loop all faces

//            const label c = patch.faceCells()[faceI];

//            // zeile 1
//            diag[c](dirI + 0, dirJ + XX) += ic[faceI].component(symmTensor::XX) * sn[faceI].x();
//            diag[c](dirI + 0, dirJ + XY) += ic[faceI].component(symmTensor::XY) * sn[faceI].y();
//            diag[c](dirI + 0, dirJ + XZ) += ic[faceI].component(symmTensor::XZ) * sn[faceI].z();

//            B[c](dirI + 0) -=   bc[faceI].component(symmTensor::XX) * sn[faceI].x()
//                              + bc[faceI].component(symmTensor::XY) * sn[faceI].y()
//                              + bc[faceI].component(symmTensor::XZ) * sn[faceI].z();

//            // zeile 2
//            diag[c](dirI + 1, dirJ + YX) += ic[faceI].component(symmTensor::XY) * sn[faceI].x();
//            diag[c](dirI + 1, dirJ + YY) += ic[faceI].component(symmTensor::YY) * sn[faceI].y();
//            diag[c](dirI + 1, dirJ + YZ) += ic[faceI].component(symmTensor::YZ) * sn[faceI].z();

//            B[c](dirI + 1) -=   bc[faceI].component(symmTensor::XY) * sn[faceI].x()
//                              + bc[faceI].component(symmTensor::YY) * sn[faceI].y()
//                              + bc[faceI].component(symmTensor::YZ) * sn[faceI].z();

//            // zeile 3
//            diag[c](dirI + 2, dirJ + ZX) += ic[faceI].component(symmTensor::XZ) * sn[faceI].x();
//            diag[c](dirI + 2, dirJ + ZY) += ic[faceI].component(symmTensor::YZ) * sn[faceI].y();
//            diag[c](dirI + 2, dirJ + ZZ) += ic[faceI].component(symmTensor::ZZ) * sn[faceI].z();

//            B[c](dirI + 2) -=   bc[faceI].component(symmTensor::XZ) * sn[faceI].x()
//                              + bc[faceI].component(symmTensor::YZ) * sn[faceI].y()
//                              + bc[faceI].component(symmTensor::ZZ) * sn[faceI].z();

//        }

//    }

}

template<class vectorType>
void gaussDiv(volVectorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("div(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    for(int i=0; i<vf.mesh().owner().size(); i++){

        int o = vf.mesh().owner()[i];
        int n = vf.mesh().neighbour()[i];
        scalar w = weights.internalField()[i];
        vector s = vf.mesh().Sf()[i];

        diag[o](dirI + 0, dirJ + 0) += f*w*s.x() * volInv[o];
        diag[o](dirI + 0, dirJ + 1) += f*w*s.y() * volInv[o];
        diag[o](dirI + 0, dirJ + 2) += f*w*s.z() * volInv[o];

        diag[n](dirI + 0, dirJ + 0) -= f*(1-w)*s.x() * volInv[n];
        diag[n](dirI + 0, dirJ + 1) -= f*(1-w)*s.y() * volInv[n];
        diag[n](dirI + 0, dirJ + 2) -= f*(1-w)*s.z() * volInv[n];

        upper[i](dirI + 0, dirJ + 0) += f*(1-w)*s.x() * volInv[o];
        upper[i](dirI + 0, dirJ + 1) += f*(1-w)*s.y() * volInv[o];
        upper[i](dirI + 0, dirJ + 2) += f*(1-w)*s.z() * volInv[o];

        lower[i](dirI + 0, dirJ + 0) -= f*w*s.x() * volInv[n];
        lower[i](dirI + 0, dirJ + 1) -= f*w*s.y() * volInv[n];
        lower[i](dirI + 0, dirJ + 2) -= f*w*s.z() * volInv[n];

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<vector>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        if (patchVolField.type() == "myImplicitVelocity"){

            //notImplemented("myImplicitVelocity für div(vectorField) noch nicht implementiert!!!");


            tmp<tensorField> ticd = patchVolField.valueInternalCoeffsTensor(weightsPatchVolField);
            tmp<tensorField> ticu = patchVolField.valueInternalCoeffsUpperTensor(weightsPatchVolField);
            tmp<tensorField> ticl = patchVolField.valueInternalCoeffsLowerTensor(weightsPatchVolField);
            const tensorField& icd = ticd();              // internal coefficient
            const tensorField& icu = ticu();              // upper coefficient
            const tensorField& icl = ticl();              // lower coefficient

            const labelList& faceCells = patchVolField.patch().faceCells();
            const vectorField& sn = patchVolField.patch().nf();
            const labelList& secondFaces = patchVolField.secondFaces();

            forAll(patchVolField, faceI){               // loop all faces

                const label c = faceCells[faceI];
                const label secondFace = secondFaces[faceI];
                const vector n = sn[faceI];
                const tensor coeffsD = icd[faceI];
                const tensor coeffsU = icu[faceI];
                const tensor coeffsL = icl[faceI];


                diag[c](dirI + 0, dirJ + 0) += f*coeffsD.xx()*n.x() + f*coeffsD.yx()*n.y() + f*coeffsD.zx()*n.z();
                diag[c](dirI + 0, dirJ + 1) += f*coeffsD.xy()*n.x() + f*coeffsD.yy()*n.y() + f*coeffsD.zy()*n.z();
                diag[c](dirI + 0, dirJ + 2) += f*coeffsD.xz()*n.x() + f*coeffsD.yz()*n.y() + f*coeffsD.zz()*n.z();

                upper[secondFace](dirI + 0, dirJ + 0) += f*coeffsU.xx()*n.x() + f*coeffsU.yx()*n.y() + f*coeffsU.zx()*n.z();
                upper[secondFace](dirI + 0, dirJ + 1) += f*coeffsU.xy()*n.x() + f*coeffsU.yy()*n.y() + f*coeffsU.zy()*n.z();
                upper[secondFace](dirI + 0, dirJ + 2) += f*coeffsU.xz()*n.x() + f*coeffsU.yz()*n.y() + f*coeffsU.zz()*n.z();

                lower[secondFace](dirI + 0, dirJ + 0) += f*coeffsL.xx()*n.x() + f*coeffsL.yx()*n.y() + f*coeffsL.zx()*n.z();
                lower[secondFace](dirI + 0, dirJ + 1) += f*coeffsL.xy()*n.x() + f*coeffsL.yy()*n.y() + f*coeffsL.zy()*n.z();
                lower[secondFace](dirI + 0, dirJ + 2) += f*coeffsL.xz()*n.x() + f*coeffsL.yz()*n.y() + f*coeffsL.zz()*n.z();

            }


            //Hier noch füllen!!!

        }
        else if(patchVolField.type() == "implicitExtrapolation"){
            notImplemented("implicitExtrapolation für div(vectorField) noch nicht implementiert!!!");
        }
        else{

            tmp<vectorField> tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
            tmp<vectorField> tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
            const vectorField& ic = tic();              // internal coefficient
            const vectorField& bc = tbc();              // boundary coefficient

            const fvPatch& patch = patchVolField.patch();   // reference to patch

            const vectorField& sn = patch.Sf();         // patch normals

            forAll(patchVolField, faceI){               // loop all faces

                const label c = patch.faceCells()[faceI];
                const vector v = cmptMultiply(ic[faceI],sn[faceI]);

                diag[c](dirI + 0, dirJ + 0) += f*v.x() * volInv[c];
                diag[c](dirI + 0, dirJ + 1) += f*v.y() * volInv[c];
                diag[c](dirI + 0, dirJ + 2) += f*v.z() * volInv[c];

                B[c](dirI) -= f*(bc[faceI] & sn[faceI]) * volInv[c];
            }

        }
    }

}

template<class vectorType>
void gaussDiv(volSymmTensorField& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("div(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    const int XX = 0;
    const int XY = 1;
    const int XZ = 2;
    const int YY = 3;
    const int YZ = 4;
    const int ZZ = 5;

    const int YX = XY;
    const int ZX = XZ;
    const int ZY = YZ;

    for(int i=0; i<vf.mesh().owner().size(); i++){

        int o = vf.mesh().owner()[i];
        int n = vf.mesh().neighbour()[i];
        scalar w = weights.internalField()[i];
        vector s = vf.mesh().Sf()[i];

        // zeile 1
        diag[o](dirI + 0, dirJ + XX) += f*w*s.x() * volInv[o];
        diag[o](dirI + 0, dirJ + XY) += f*w*s.y() * volInv[o];
        diag[o](dirI + 0, dirJ + XZ) += f*w*s.z() * volInv[o];

        diag[n](dirI + 0, dirJ + XX) -= f*(1-w)*s.x() * volInv[n];
        diag[n](dirI + 0, dirJ + XY) -= f*(1-w)*s.y() * volInv[n];
        diag[n](dirI + 0, dirJ + XZ) -= f*(1-w)*s.z() * volInv[n];

        upper[i](dirI + 0, dirJ + XX) += f*(1-w)*s.x() * volInv[o];
        upper[i](dirI + 0, dirJ + XY) += f*(1-w)*s.y() * volInv[o];
        upper[i](dirI + 0, dirJ + XZ) += f*(1-w)*s.z() * volInv[o];

        lower[i](dirI + 0, dirJ + XX) -= f*w*s.x() * volInv[n];
        lower[i](dirI + 0, dirJ + XY) -= f*w*s.y() * volInv[n];
        lower[i](dirI + 0, dirJ + XZ) -= f*w*s.z() * volInv[n];

        // zeile 2
        diag[o](dirI + 1, dirJ + YX) += f*w*s.x() * volInv[o];
        diag[o](dirI + 1, dirJ + YY) += f*w*s.y() * volInv[o];
        diag[o](dirI + 1, dirJ + YZ) += f*w*s.z() * volInv[o];

        diag[n](dirI + 1, dirJ + YX) -= f*(1-w)*s.x() * volInv[n];
        diag[n](dirI + 1, dirJ + YY) -= f*(1-w)*s.y() * volInv[n];
        diag[n](dirI + 1, dirJ + YZ) -= f*(1-w)*s.z() * volInv[n];

        upper[i](dirI + 1, dirJ + YX) += f*(1-w)*s.x() * volInv[o];
        upper[i](dirI + 1, dirJ + YY) += f*(1-w)*s.y() * volInv[o];
        upper[i](dirI + 1, dirJ + YZ) += f*(1-w)*s.z() * volInv[o];

        lower[i](dirI + 1, dirJ + YX) -= f*w*s.x() * volInv[n];
        lower[i](dirI + 1, dirJ + YY) -= f*w*s.y() * volInv[n];
        lower[i](dirI + 1, dirJ + YZ) -= f*w*s.z() * volInv[n];

        // zeile 3
        diag[o](dirI + 2, dirJ + ZX) += f*w*s.x() * volInv[o];
        diag[o](dirI + 2, dirJ + ZY) += f*w*s.y() * volInv[o];
        diag[o](dirI + 2, dirJ + ZZ) += f*w*s.z() * volInv[o];

        diag[n](dirI + 2, dirJ + ZX) -= f*(1-w)*s.x() * volInv[n];
        diag[n](dirI + 2, dirJ + ZY) -= f*(1-w)*s.y() * volInv[n];
        diag[n](dirI + 2, dirJ + ZZ) -= f*(1-w)*s.z() * volInv[n];

        upper[i](dirI + 2, dirJ + ZX) += f*(1-w)*s.x() * volInv[o];
        upper[i](dirI + 2, dirJ + ZY) += f*(1-w)*s.y() * volInv[o];
        upper[i](dirI + 2, dirJ + ZZ) += f*(1-w)*s.z() * volInv[o];

        lower[i](dirI + 2, dirJ + ZX) -= f*w*s.x() * volInv[n];
        lower[i](dirI + 2, dirJ + ZY) -= f*w*s.y() * volInv[n];
        lower[i](dirI + 2, dirJ + ZZ) -= f*w*s.z() * volInv[n];
    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<symmTensor>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<symmTensorField> tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<symmTensorField> tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const symmTensorField& ic = tic();              // internal coefficient
        const symmTensorField& bc = tbc();              // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        const vectorField& sn = patch.Sf();         // patch normals

        forAll(patchVolField, faceI){               // loop all faces

            const label c = patch.faceCells()[faceI];

            // zeile 1
            diag[c](dirI + 0, dirJ + XX) += f*ic[faceI].component(symmTensor::XX) * sn[faceI].x() * volInv[c];
            diag[c](dirI + 0, dirJ + XY) += f*ic[faceI].component(symmTensor::XY) * sn[faceI].y() * volInv[c];
            diag[c](dirI + 0, dirJ + XZ) += f*ic[faceI].component(symmTensor::XZ) * sn[faceI].z() * volInv[c];

            B[c](dirI + 0) -=   f*bc[faceI].component(symmTensor::XX) * sn[faceI].x() * volInv[c]
                              + f*bc[faceI].component(symmTensor::XY) * sn[faceI].y() * volInv[c]
                              + f*bc[faceI].component(symmTensor::XZ) * sn[faceI].z() * volInv[c];

            // zeile 2
            diag[c](dirI + 1, dirJ + YX) += f*ic[faceI].component(symmTensor::XY) * sn[faceI].x() * volInv[c];
            diag[c](dirI + 1, dirJ + YY) += f*ic[faceI].component(symmTensor::YY) * sn[faceI].y() * volInv[c];
            diag[c](dirI + 1, dirJ + YZ) += f*ic[faceI].component(symmTensor::YZ) * sn[faceI].z() * volInv[c];

            B[c](dirI + 1) -=   bc[faceI].component(symmTensor::XY) * sn[faceI].x() * volInv[c]
                              + bc[faceI].component(symmTensor::YY) * sn[faceI].y() * volInv[c]
                              + bc[faceI].component(symmTensor::YZ) * sn[faceI].z() * volInv[c];

            // zeile 3
            diag[c](dirI + 2, dirJ + ZX) += f*ic[faceI].component(symmTensor::XZ) * sn[faceI].x() * volInv[c];
            diag[c](dirI + 2, dirJ + ZY) += f*ic[faceI].component(symmTensor::YZ) * sn[faceI].y() * volInv[c];
            diag[c](dirI + 2, dirJ + ZZ) += f*ic[faceI].component(symmTensor::ZZ) * sn[faceI].z() * volInv[c];

            B[c](dirI + 2) -=   f*bc[faceI].component(symmTensor::XZ) * sn[faceI].x() * volInv[c]
                              + f*bc[faceI].component(symmTensor::YZ) * sn[faceI].y() * volInv[c]
                              + f*bc[faceI].component(symmTensor::ZZ) * sn[faceI].z() * volInv[c];

        }

    }

}

}
