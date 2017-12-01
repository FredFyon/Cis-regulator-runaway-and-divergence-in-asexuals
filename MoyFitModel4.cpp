#include "MoyFitModel4.h"
#include "SelectionModel4.h"
#include <string>
#include <vector>

using namespace std;

double MoyFitModel4(double Ei, vector<double> PopEChrom1, vector<double> PopSChrom1, vector<double> PopEChrom2, vector<double> PopSChrom2, double h, int Npop, string FIT, double I)
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
            WMOY += SelectionModel4(Ei,PopEChrom1[indFit],PopSChrom1[indFit],PopEChrom2[indFit],PopSChrom2[indFit],h,I).getx1();
        }
    }

    WMOY /= Npop;

    return WMOY;
}
