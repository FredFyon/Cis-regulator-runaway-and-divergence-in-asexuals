#include "SelectionModel6.h"
#include <string>
#include <cmath>
#include <iostream>
#include "math.h"

using namespace std;

double SelectionModel6(double Ei, double DChrom1, double EChrom1, double SChrom1, double DChrom2, double EChrom2, double SChrom2, double h, double I, string Modifier)
{
    double S = max(SChrom1,SChrom2);

    double y,Wes;
    double E;

    if(Modifier == "QProt")
    {
        E = DChrom1*DChrom2*(log10(EChrom1)+log10(EChrom2));
    }
    else if(Modifier == "ExcesQProt")
    {
        E = 2*log10(Ei) + (log10(EChrom1)+log10(EChrom2)-2*log10(Ei))*(DChrom1*DChrom2);
    }

    if((EChrom1 == 0) && (EChrom2 == 0))
    {
        cout << "Promoters attractivity for transcription factors has collapsed. Model can't calculate individual fitness." << endl;
    }
    else
    {
        if(S == SChrom1)
        {
            y = EChrom1/(EChrom1+EChrom2);
            Wes = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
        }
        else
        {
            y = EChrom2/(EChrom1+EChrom2);
            Wes = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
        }
    }

    double Wconstraint;

    Wconstraint = exp(-I*pow(E-2*log10(Ei),2));

    double W = Wes*Wconstraint;

    return W;
}
