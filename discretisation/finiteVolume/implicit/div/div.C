#include"fvCFD.H"
#include "../useful/interpolationWeights.H"

namespace myFvm{


template<class Type, class vectorType>
void gaussDiv(GeometricField<Type, fvPatchField, volMesh>& vf, const scalar f, BlockLduMatrix<vectorType>& M, Field<vectorType>& B, const int dirI, const int dirJ, const bool updateOnlyRHS = false){

    typename CoeffField<vectorType>::squareTypeField& diag = M.diag().asSquare();
    typename CoeffField<vectorType>::squareTypeField& upper = M.upper().asSquare();
    typename CoeffField<vectorType>::squareTypeField& lower = M.lower().asSquare();

    const word name("div(" + vf.name() + ")");
    tmp<surfaceScalarField> tweights = interpolWeights(vf, name);
    const surfaceScalarField& weights = tweights();

    tmp<scalarField> tvol = vf.mesh().V();
    const scalarField volInv = 1/tvol();

    int dimBild;
    int components[3^(Type::rank)/3][3];

    const string fieldType = Type::typeName;
    if (fieldType == "vector"){

        dimBild = 1;
        components[0][0] = 0;       // X
        components[0][1] = 1;       // Y
        components[0][2] = 2;       // Z

    }
    else if(fieldType == "symmTensor"){

        dimBild = 3;
        components[0][0] = 0;       // XX
        components[0][1] = 1;       // XY
        components[0][2] = 2;       // XZ
        components[1][0] = 1;       // YX
        components[1][1] = 3;       // YY
        components[1][2] = 4;       // YZ
        components[2][0] = 2;       // ZX
        components[2][1] = 4;       // ZY
        components[2][2] = 5;       // ZZ

    }
    else{

        Info << "notImplemented: myFvm::gaussDiv(" << Type::typeName << ")" << endl;
        notImplemented("siehe Infomeldung")

    }


    if(!updateOnlyRHS){

        for(int i=0; i<vf.mesh().owner().size(); i++){

            int o = vf.mesh().owner()[i];
            int n = vf.mesh().neighbour()[i];
            scalar w = weights.internalField()[i];
            vector s = vf.mesh().Sf()[i];

            for(int j=0; j<dimBild; j++){

                diag[o](dirI + j, dirJ + components[j][0]) += f*w*s.x() * volInv[o];
                diag[o](dirI + j, dirJ + components[j][1]) += f*w*s.y() * volInv[o];
                diag[o](dirI + j, dirJ + components[j][2]) += f*w*s.z() * volInv[o];

                diag[n](dirI + j, dirJ + components[j][0]) -= f*(1-w)*s.x() * volInv[n];
                diag[n](dirI + j, dirJ + components[j][1]) -= f*(1-w)*s.y() * volInv[n];
                diag[n](dirI + j, dirJ + components[j][2]) -= f*(1-w)*s.z() * volInv[n];

                upper[i](dirI + j, dirJ + components[j][0]) += f*(1-w)*s.x() * volInv[o];
                upper[i](dirI + j, dirJ + components[j][1]) += f*(1-w)*s.y() * volInv[o];
                upper[i](dirI + j, dirJ + components[j][2]) += f*(1-w)*s.z() * volInv[o];

                lower[i](dirI + j, dirJ + components[j][0]) -= f*w*s.x() * volInv[n];
                lower[i](dirI + j, dirJ + components[j][1]) -= f*w*s.y() * volInv[n];
                lower[i](dirI + j, dirJ + components[j][2]) -= f*w*s.z() * volInv[n];

            }

        }

    }

    // CONTRIBUTION BOUNDARY FIELD ...
    //
    //      vfboundFace = ic * vfboundCell  +  bc

    vf.boundaryField().updateCoeffs();
    forAll(vf.boundaryField(), patchI){ // loop all patches

        const fvPatchField<Type>& patchVolField = vf.boundaryField()[patchI];
        const fvsPatchScalarField& weightsPatchVolField = weights.boundaryField()[patchI];  // weight that coeffs get multiplied with

        tmp<Field<Type> > tic = patchVolField.valueInternalCoeffs(weightsPatchVolField);
        tmp<Field<Type> > tbc = patchVolField.valueBoundaryCoeffs(weightsPatchVolField);
        const Field<Type>& ic = tic();              // internal coefficient
        const Field<Type>& bc = tbc();              // boundary coefficient

        const fvPatch& patch = patchVolField.patch();   // reference to patch

        const vectorField& sn = patch.Sf();         // patch normals

        forAll(patchVolField, faceI){               // loop all faces

            const label c = patch.faceCells()[faceI];

            if(!updateOnlyRHS){

                for(int j=0; j<dimBild; j++){

                    diag[c](dirI + 0, dirJ + components[j][0]) += f*ic[faceI].component(components[j][0]) * sn[faceI].x() * volInv[c];
                    diag[c](dirI + 0, dirJ + components[j][1]) += f*ic[faceI].component(components[j][1]) * sn[faceI].y() * volInv[c];
                    diag[c](dirI + 0, dirJ + components[j][2]) += f*ic[faceI].component(components[j][2]) * sn[faceI].z() * volInv[c];

                }

            }

            for(int j=0; j<dimBild; j++){

                B[c](dirI + j) -=   f*bc[faceI].component(components[j][0]) * sn[faceI].x() * volInv[c]
                                  + f*bc[faceI].component(components[j][1]) * sn[faceI].y() * volInv[c]
                                  + f*bc[faceI].component(components[j][2]) * sn[faceI].z() * volInv[c];

            }

        }

    }

}

}
