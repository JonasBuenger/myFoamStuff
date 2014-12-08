
#include"laplacian.H"
#include"fvCFD.H"
#include"fvMesh.H"

namespace myFvm{

subMatrix<subBlock<1,1>, subBlock<1,1> > laplacian(volScalarField& vf){

    subMatrix<subBlock<1,1>, subBlock<1,1> > matrix(vf.mesh());

    Foam::fvScalarMatrix m = fvm::laplacian(vf);
    m.addBoundaryDiag(m.diag(),0);
    m.addBoundarySource(m.source(),false);

    forAll(m.diag(), i){
        matrix.diag(i).element(0,0) = m.diag()[i];
    }

    forAll(m.source(), i){
        matrix.source(i).element(0,0) = m.source()[i];
    }

    // prüfen, ob upper überhaupt notwendig
    forAll(m.upper(), i){
        matrix.upper(i).element(0,0) = m.upper()[i];
    }

    // prüfen, ob lower überhaupt notwendig
    forAll(m.lower(), i){
        matrix.lower(i).element(0,0) = m.lower()[i];
    }

    return matrix;

}

subMatrix<subBlock<3,3>, subBlock<3,1> > laplacian(volVectorField &vf){

    subMatrix<subBlock<3,3>, subBlock<3,1> > matrix(vf.mesh());

    Foam::fvVectorMatrix m = fvm::laplacian(vf);
    m.addBoundaryDiag(m.diag(),0);
    m.addBoundarySource(m.source(),false);

    forAll(m.diag(), i){
        matrix.diag(i).element(0,0) = m.diag()[i];
        matrix.diag(i).element(1,1) = m.diag()[i];
        matrix.diag(i).element(2,2) = m.diag()[i];
    }

    forAll(m.source(), i){
        matrix.source(i).element(0,0) = m.source()[i].component(0);
        matrix.source(i).element(1,0) = m.source()[i].component(1);
        matrix.source(i).element(2,0) = m.source()[i].component(2);
    }

    // prüfen, ob upper überhaupt notwendig
    forAll(m.upper(), i){
        matrix.upper(i).element(0,0) = m.upper()[i];
        matrix.upper(i).element(1,1) = m.upper()[i];
        matrix.upper(i).element(2,2) = m.upper()[i];
    }

    // prüfen, ob lower überhaupt notwendig
    forAll(m.lower(), i){
        matrix.lower(i).element(0,0) = m.lower()[i];
        matrix.lower(i).element(1,1) = m.lower()[i];
        matrix.lower(i).element(2,2) = m.lower()[i];
    }

    return matrix;

}

}
