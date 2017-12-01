#include "SelectionModel7.h"
#include <iostream>
#include <math.h>

using namespace std;

double SelectionModel7(double Ei, double E1Chrom1, double E2Chrom1, double SChrom1, double E1Chrom2, double E2Chrom2, double SChrom2, double h, double I)
{
    double S = max(SChrom1,SChrom2);

    double y,WA;

    double EChrom1 = E1Chrom1+E2Chrom1;
    double EChrom2 = E1Chrom2+E2Chrom2;

    if(S == SChrom1)
    {
        y = EChrom1/(EChrom1+EChrom2);
        WA = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
    }
    else
    {
        y = EChrom2/(EChrom1+EChrom2);
        WA = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
    }

    double E = log10(E1Chrom1)+log10(E2Chrom1)+log10(E1Chrom2)+log10(E2Chrom2);
    double Wconstraint;

    Wconstraint = exp(-I*pow(E-4*log10(Ei),2));
    //cout << WA << " "  << Wconstraint << endl;

    double W = WA*Wconstraint;

    //cout << WA << " " << Wconstraint << endl;

    return W;
}

