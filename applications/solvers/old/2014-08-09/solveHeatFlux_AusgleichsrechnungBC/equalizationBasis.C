#include "equalizationBasis.h"
#include "volFields.H"


Foam::scalar Foam::equalizationBasis::c(const vector center, const vector point){
    return 1;
}

Foam::scalar Foam::equalizationBasis::x(const vector center, const vector point){
    return (point-center).x();
}

Foam::scalar Foam::equalizationBasis::y(const vector center, const vector point){
    return (point-center).y();
}

Foam::scalar Foam::equalizationBasis::z(const vector center, const vector point){
    return (point-center).z();
}

Foam::equalizationBasis::equalizationBasis(const int t)
{
    switch(t){
    case 0:
        _type = "constant";
        _dim = 1;
        break;
    case 1:
        _type = "linear";
        _dim = 4;
        break;
    case 2:
        _type = "bilinear";
        _dim = 7;
        break;
    case 3:
        _type = "quadratic";
        _dim = 10;
        break;
    default:
        Info << "Ungültiger Polynomtyp" << endl;
        break;
    }
}

Foam::scalar Foam::equalizationBasis::eval(const vector center, const vector point, const int numBasisFcn){

    // FEHLERABFRAGE: passt numBasisFcn zum type?

    scalar res = 0;

    switch(numBasisFcn){
    case 0:
        res = c(center, point);
        break;
    case 1:
        res = x(center, point);
        break;
    case 2:
        res = y(center, point);
        break;
    case 3:
        res = z(center, point);
        break;
    case 4:
        res = x(center, point)*y(center,point);
        break;
    case 5:
        res = x(center, point)*z(center,point);
        break;
    case 6:
        res = y(center, point)*z(center,point);
        break;
    }

    if (res < 1e-10 && res > -1e-10)
        res = 0;

    return res;

}

Foam::scalar Foam::equalizationBasis::dim(){
    return _dim;
}
