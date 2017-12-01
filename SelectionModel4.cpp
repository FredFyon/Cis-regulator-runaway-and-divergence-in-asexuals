#include "SelectionModel4.h"
#include "SelectionModel1.h"
#include <string>
#include <cmath>
#include "Vecteur.h"

using namespace std;

Vecteur SelectionModel4(double Ei, double EChrom1, double SChrom1, double EChrom2, double SChrom2, double h, double I)
{
    double Wes = SelectionModel1(EChrom1,SChrom1,EChrom2,SChrom2,h).getx1();
    double H = SelectionModel1(EChrom1,SChrom1,EChrom2,SChrom2,h).getx2();

    double E = log10(EChrom1) + log10(EChrom2);
    double Wconstraint;

    Wconstraint = exp(-I*pow(E-2*log10(Ei),2));

    double W = Wes*Wconstraint;

    Vecteur Res(W,H,0.,0.,0.,0.,0.,0.,0.,0.);

    return Res;
}
