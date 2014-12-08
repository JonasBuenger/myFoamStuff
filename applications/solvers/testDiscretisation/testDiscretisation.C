/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    blockCoupledScalarTransportFoam

Description
    Solves two coupled transport equations in a block-coupled manner

        1) transport equation for a passive scalar
        2) diffusion only

    This resembles heat exchanging flow through a porous medium

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockLduSolvers.H"
#include "blockVectorNMatrices.H"
#include "blockMatrixTools.H"
#include "BlockSolverPerformance.H"

#include "../../../discretisation/finiteVolume/myFvm.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nTESTING OPERATORS \n" << endl;

    while(runTime.loop()) //for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool Sp = true;
        bool div = true;
        bool grad = true;
        bool laplacian = true;

        // initialize fields
        forAll(sc, i){

            //sc[i] = i+1;
            sc[i] = drand48();
            s[i] = sc[i];
            vc[i] = vector(drand48(), drand48(), drand48());
            //vc[i] = vector(i+1, (i+1)*10, (i+1)*100);
            v[i] = vc[i];
            sTc[i] = symmTensor(drand48(), drand48(), drand48(), drand48(), drand48(), drand48());
            //sTc[i] = symmTensor(i+1, (i+1)*10, (i+1)*100, (i+1)*1000, (i+1)*10000, (i+1)*100000);
            sT[i] = sTc[i];

        }

        // Test weather all fields are equal
        if (max(sc.internalField()-s.internalField()) > 0.0)
            Info << "sc-s != 0" << endl;

        if (max(vc.internalField().component(0)-v.internalField().component(0)) > 0.0)
            Info << "vc.x - v.x != 0" << endl;
        if (max(vc.internalField().component(1)-v.internalField().component(1)) > 0.0)
            Info << "vc.z - v.z != 0" << endl;
        if (max(vc.internalField().component(2)-v.internalField().component(2)) > 0.0)
            Info << "vc.z - v.z != 0" << endl;

        if (max(sTc.internalField().component(0)-sT.internalField().component(0)) > 0.0)
            Info << "sTc.xx - sT.xx != 0" << endl;
        if (max(sTc.internalField().component(1)-sT.internalField().component(1)) > 0.0)
            Info << "sTc.xy - sT.xy != 0" << endl;
        if (max(sTc.internalField().component(2)-sT.internalField().component(2)) > 0.0)
            Info << "sTc.xz - sT.xz != 0" << endl;
        if (max(sTc.internalField().component(3)-sT.internalField().component(3)) > 0.0)
            Info << "sTc.yy - sT.yy != 0" << endl;
        if (max(sTc.internalField().component(4)-sT.internalField().component(4)) > 0.0)
            Info << "sTc.yz - sT.yz != 0" << endl;
        if (max(sTc.internalField().component(5)-sT.internalField().component(5)) > 0.0)
            Info << "sTc.zz - sT.zz != 0" << endl;

        sc.correctBoundaryConditions();
        vc.correctBoundaryConditions();
        sTc.correctBoundaryConditions();

            // vector1
//            BlockLduMatrix<vector1> M1(mesh);
//            Field<vector1> B1(mesh.nCells(), vector1::zero);
            // vector3
            BlockLduMatrix<vector3> M3(mesh);
            Field<vector3> B3(mesh.nCells(), vector3::zero);
            Field<vector3> MB3(mesh.nCells(), vector3::zero);
            // vector6
            BlockLduMatrix<vector6> M6(mesh);
            Field<vector6> B6(mesh.nCells(), vector6::zero);
            Field<vector6> MB6(mesh.nCells(), vector6::zero);

            Field<vector3> Ax3(X3.size());
            Field<vector6> Ax6(X6.size());


            volScalarField f = s;
            forAll(f, i){
                f[i] = 1;
            }


            // ----------------GRADIENT---------------//
            if (grad){

                Info << endl;

                // grad(scalar)
                M3.diag().asSquare() = tensor3::zero;
                M3.upper().asSquare() = tensor3::zero;
                M3.lower().asSquare() = tensor3::zero;
                B3 = vector3::zero;
                MB3 = vector3::zero;
                s = sc;

                blockMatrixTools::blockInsert(0,sc.component(0),X3);

                tmp<volVectorField> gradS = fvc::grad(sc);
                blockMatrixTools::blockInsert(0,gradS().component(0),B3);
                blockMatrixTools::blockInsert(1,gradS().component(1),B3);
                blockMatrixTools::blockInsert(2,gradS().component(2),B3);

                myFvm::gaussGrad(s,1,M3,MB3,0,0);

                M3.Amul(Ax3, X3.internalField());

                Info << "Max-Error in gaussGrad(scalar): " << max(mag((Ax3-MB3)-B3)) << endl;

                // grad(vector)
                M6.diag().asSquare() = tensor6::zero;
                M6.upper().asSquare() = tensor6::zero;
                M6.lower().asSquare() = tensor6::zero;
                B6 = vector6::zero;
                MB6 = vector6::zero;
                v = vc;

                const label XX = 0;
                const label XY = 1;
                const label XZ = 2;
                const label YX = 3;
                const label YY = 4;
                const label YZ = 5;
                const label ZX = 6;
                const label ZY = 7;
                const label ZZ = 8;

                blockMatrixTools::blockInsert(0,vc.component(0),X6);
                blockMatrixTools::blockInsert(1,vc.component(1),X6);
                blockMatrixTools::blockInsert(2,vc.component(2),X6);

                tmp<volTensorField> gradV = fvc::grad(vc);
                blockMatrixTools::blockInsert(0,0.5*(gradV().component(XX)+gradV().component(XX)),B6);
                blockMatrixTools::blockInsert(1,0.5*(gradV().component(XY)+gradV().component(YX)),B6);
                blockMatrixTools::blockInsert(2,0.5*(gradV().component(XZ)+gradV().component(ZX)),B6);
                blockMatrixTools::blockInsert(3,0.5*(gradV().component(YY)+gradV().component(YY)),B6);
                blockMatrixTools::blockInsert(4,0.5*(gradV().component(YZ)+gradV().component(ZY)),B6);
                blockMatrixTools::blockInsert(5,0.5*(gradV().component(ZZ)+gradV().component(ZZ)),B6);

                myFvm::gaussSymmGrad(v,1,M6,MB6,0,0);

                M6.Amul(Ax6,X6.internalField());

                Info << "Max-Error in gaussSymmGrad(vector): " << max(mag((Ax6-MB6)-B6)) << endl;

            }

            // ----------------DIVERGENZ---------------//
            if (div){

                Info << endl;

                // vector
                M3.diag().asSquare() = tensor3::zero;
                M3.upper().asSquare() = tensor3::zero;
                M3.lower().asSquare() = tensor3::zero;
                B3 = vector3::zero;
                MB3 = vector3::zero;
                v = vc;

                tmp<volScalarField> divV = fvc::div(vc);
                blockMatrixTools::blockInsert(0,divV().component(0),B3);

                myFvm::gaussDiv(v,1,M3,B3,0,0);

                blockMatrixTools::blockInsert(0,vc.internalField().component(0),X3);
                blockMatrixTools::blockInsert(1,vc.internalField().component(1),X3);
                blockMatrixTools::blockInsert(2,vc.internalField().component(2),X3);

                M3.Amul(Ax3,X3.internalField());

                Info << "Max-Error in gaussDiv(vector): "
                     << max(mag((Ax3.component(0)-MB3.component(0))-B3.component(0))) << endl;


                // symmTensor
                M6.diag().asSquare() = tensor6::zero;
                M6.upper().asSquare() = tensor6::zero;
                M6.lower().asSquare() = tensor6::zero;
                B6 = vector6::zero;
                MB6 = vector6::zero;
                sT = sTc;

                tmp<volVectorField> divsT = fvc::div(sTc);
                blockMatrixTools::blockInsert(0,divsT().component(0),B6);
                blockMatrixTools::blockInsert(1,divsT().component(1),B6);
                blockMatrixTools::blockInsert(2,divsT().component(2),B6);

                blockMatrixTools::blockInsert(0,sTc.internalField().component(0),X6);
                blockMatrixTools::blockInsert(1,sTc.internalField().component(1),X6);
                blockMatrixTools::blockInsert(2,sTc.internalField().component(2),X6);
                blockMatrixTools::blockInsert(3,sTc.internalField().component(3),X6);
                blockMatrixTools::blockInsert(4,sTc.internalField().component(4),X6);
                blockMatrixTools::blockInsert(5,sTc.internalField().component(5),X6);

                myFvm::gaussDiv(sT,1,M6,MB6,0,0);

                M6.Amul(Ax6,X6.internalField());

                Info << "Max-Error in gaussDiv(symmTensor): " << mag(vector(
                     max((Ax6.component(0)-MB6.component(0))-B6.component(0)),
                     max((Ax6.component(1)-MB6.component(1))-B6.component(1)),
                     max((Ax6.component(2)-MB6.component(2))-B6.component(2)) )) << endl;
            }

            // ----------------SP---------------//
            if (Sp){

                Info << endl;

                // vector
                M3.diag().asSquare() = tensor3::zero;
                M3.upper().asSquare() = tensor3::zero;
                M3.lower().asSquare() = tensor3::zero;
                B3 = vector3::zero;
                MB3 = vector3::zero;
                v = vc;

                blockMatrixTools::blockInsert(0,v.internalField().component(0),X3);
                blockMatrixTools::blockInsert(1,v.internalField().component(1),X3);
                blockMatrixTools::blockInsert(2,v.internalField().component(2),X3);

                blockMatrixTools::blockInsert(0,sT.internalField().component(0),X6);
                blockMatrixTools::blockInsert(1,sT.internalField().component(1),X6);
                blockMatrixTools::blockInsert(2,sT.internalField().component(2),X6);
                blockMatrixTools::blockInsert(3,sT.internalField().component(3),X6);
                blockMatrixTools::blockInsert(4,sT.internalField().component(4),X6);
                blockMatrixTools::blockInsert(5,sT.internalField().component(5),X6);

                tmp<volVectorField> SpV = fvc::Sp(f ,vc);
                blockMatrixTools::blockInsert(0,SpV().component(0),B3);
                blockMatrixTools::blockInsert(1,SpV().component(1),B3);
                blockMatrixTools::blockInsert(2,SpV().component(2),B3);

                myFvm::Sp(v,1,M3,MB3,0,0);

                M3.Amul(Ax3,X3.internalField());

                Info << "Max-Error in Sp(vector): " << max(mag((Ax3-MB3)-B3)) << endl;

                // symmTensor
                M6.diag().asSquare() = tensor6::zero;
                M6.upper().asSquare() = tensor6::zero;
                M6.lower().asSquare() = tensor6::zero;
                B6 = vector6::zero;
                MB6 = vector6::zero;
                sT = sTc;

                tmp<volSymmTensorField> SpsT = fvc::Sp(f,sTc);
                blockMatrixTools::blockInsert(0,SpsT().component(0),B6);
                blockMatrixTools::blockInsert(1,SpsT().component(1),B6);
                blockMatrixTools::blockInsert(2,SpsT().component(2),B6);
                blockMatrixTools::blockInsert(3,SpsT().component(3),B6);
                blockMatrixTools::blockInsert(4,SpsT().component(4),B6);
                blockMatrixTools::blockInsert(5,SpsT().component(5),B6);

                myFvm::Sp(sT,1,M6,MB6,0,0);

                M6.Amul(Ax6,X6.internalField());

                Info << "Max-Error in Sp(symmTensor): " << max(mag((Ax6-MB6)-B6)) << endl;

            }

            //-------------LAPLACIAN-----------//
            if (laplacian){

                Info << endl;

                // vector
                M3.diag().asSquare() = tensor3::zero;
                M3.upper().asSquare() = tensor3::zero;
                M3.lower().asSquare() = tensor3::zero;
                B3 = vector3::zero;
                MB3 = vector3::zero;
                v = vc;

                tmp<volVectorField> laplacianV = fvc::laplacian(vc);
                blockMatrixTools::blockInsert(0,laplacianV().component(0),B3);
                blockMatrixTools::blockInsert(1,laplacianV().component(1),B3);
                blockMatrixTools::blockInsert(2,laplacianV().component(2),B3);

                myFvm::gaussLaplacian(v,1,M3,MB3,0,0);

                blockMatrixTools::blockInsert(0,v.internalField().component(0),X3);
                blockMatrixTools::blockInsert(1,v.internalField().component(1),X3);
                blockMatrixTools::blockInsert(2,v.internalField().component(2),X3);

                M3.Amul(Ax3,X3.internalField());

                Info << "Max-Error in laplacian(vector): " << max(mag((Ax3-MB3)-B3)) << endl;

                // symmTensor
                M6.diag().asSquare() = tensor6::zero;
                M6.upper().asSquare() = tensor6::zero;
                M6.lower().asSquare() = tensor6::zero;
                B6 = vector6::zero;
                MB6 = vector6::zero;
                sT = sTc;

                tmp<volSymmTensorField> laplaciansT = fvc::laplacian(sTc);
                blockMatrixTools::blockInsert(0,laplaciansT().component(0),B6);
                blockMatrixTools::blockInsert(1,laplaciansT().component(1),B6);
                blockMatrixTools::blockInsert(2,laplaciansT().component(2),B6);
                blockMatrixTools::blockInsert(3,laplaciansT().component(3),B6);
                blockMatrixTools::blockInsert(4,laplaciansT().component(4),B6);
                blockMatrixTools::blockInsert(5,laplaciansT().component(5),B6);

                myFvm::gaussLaplacian(sT,1,M6,MB6,0,0);

                blockMatrixTools::blockInsert(0,sT.internalField().component(0),X6);
                blockMatrixTools::blockInsert(1,sT.internalField().component(1),X6);
                blockMatrixTools::blockInsert(2,sT.internalField().component(2),X6);
                blockMatrixTools::blockInsert(3,sT.internalField().component(3),X6);
                blockMatrixTools::blockInsert(4,sT.internalField().component(4),X6);
                blockMatrixTools::blockInsert(5,sT.internalField().component(5),X6);

                M6.Amul(Ax6,X6.internalField());

                Info << "Max-Error in laplacian(symmTensor): " << max(mag((Ax6-MB6)-B6)) << endl;

            }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
