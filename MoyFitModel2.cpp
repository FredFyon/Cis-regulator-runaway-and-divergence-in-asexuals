#include "MoyFitModel2.h"
#include "SelectionModel2.h"
#include <vector>
#include <string>

using namespace std;

double MoyFitModel2(vector<double> EChrom1a, vector<double> SChrom1a, vector<double> EChrom2a, vector<double> SChrom2a, vector<double> EChrom1b, vector<double> SChrom1b, vector<double> EChrom2b, vector<double> SChrom2b, double h, int Npop, double I, double Concaveness, string FIT, string SHAPE)
{
    double WMOY(0);
    int indFit(0);

    for(indFit = 0 ; indFit < Npop ; ++indFit)
    {
        if(FIT == "Derive")
        {
            WMOY += 1.;
        }
        else if(FIT == "Selection")
        {
            WMOY += SelectionModel2(EChrom1a[indFit],SChrom1a[indFit],EChrom2a[indFit],SChrom2a[indFit],EChrom1b[indFit],SChrom1b[indFit],EChrom2b[indFit],SChrom2b[indFit],h,I,Concaveness,SHAPE);
        }
    }

    WMOY /= Npop;

    return WMOY;
}
