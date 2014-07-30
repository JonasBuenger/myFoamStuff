#ifndef EQUALIZATIONBASIS_H
#define EQUALIZATIONBASIS_H

#include "word.H"
#include "scalar.H"
#include "vector.H"

namespace Foam{

class equalizationBasis
{
private:
    word _type;
    scalar _dim;
    scalar c(const vector center, const vector point);
    scalar x(const vector center, const vector point);
    scalar y(const vector center, const vector point);
    scalar z(const vector center, const vector point);
public:
    equalizationBasis(const int t);
    scalar dim();
    word type();
    scalar eval(const vector center, const vector point, const int numBasisFcn);
};

}
#endif // EQUALIZATIONBASIS_H
