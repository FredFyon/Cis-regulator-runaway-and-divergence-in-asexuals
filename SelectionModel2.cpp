#include "SelectionModel1.h"
#include "SelectionModel2.h"
#include <string>
#include <cmath>

using namespace std;

double SelectionModel2(double EChrom1a, double SChrom1a, double EChrom2a, double SChrom2a, double EChrom1b, double SChrom1b, double EChrom2b, double SChrom2b, double h, double I, double Concaveness, string SHAPE)
{
    double Wa = SelectionModel1(EChrom1a,SChrom1a,EChrom2a,SChrom2a,h).getx1();
    double Wb = SelectionModel1(EChrom1b,SChrom1b,EChrom2b,SChrom2b,h).getx1();

    double LogEa = log10(EChrom1a) + log10(EChrom2a);
    double LogEb = log10(EChrom1b) + log10(EChrom2b);
    double Wconstraint;

    Wconstraint = exp(-I*pow(LogEa-LogEb,2));

    double W = Wa*Wb*Wconstraint;

    return W;
}
