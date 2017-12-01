#include "SelectionModel3.h"
#include "SelectionModel1.h"
#include <string>
#include <cmath>
#include <iostream>

using namespace std;

double SelectionModel3(double TChrom1, double EChrom1, double SChrom1, double TChrom2, double EChrom2, double SChrom2, double h, double I, double Concaveness, string SHAPE)
{

    double Wes = SelectionModel1(EChrom1,SChrom1,EChrom2,SChrom2,h).getx1();

    double LogT = log10(TChrom1) + log10(TChrom2);
    double LogE = log10(EChrom1) + log10(EChrom2);
    double Wconstraint;

    Wconstraint = exp(-I*pow(LogE+LogT,2));

    double W = Wes*Wconstraint;

    //cout << Wconstraint << endl;

    return W;
}
