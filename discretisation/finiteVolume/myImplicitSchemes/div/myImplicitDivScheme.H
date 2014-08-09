#ifndef MYIMPLICITDIVSCHEME_H
#define MYIMPLICITDIVSCHEME_H

#include "fvCFD.H"

namespace Foam{

class myImplicitDivScheme{

private:
    volVectorField& vf_;
    const fvMesh& mesh_;
    bool upToDate_;
    vectorField diag_;
    vectorField upper_;
    vectorField lower_;
    scalarField source_;
    tmp<surfaceScalarField> getWeights();

public:
    myImplicitDivScheme(volVectorField&);
    void updateCoeffs();
    const vectorField& coeffsDiag();
    const vectorField& coeffsUpper();
    const vectorField& coeffsLower();
    const scalarField& coeffsSource();

};


}
#endif // MYIMPLICITDIVSCHEME_H
