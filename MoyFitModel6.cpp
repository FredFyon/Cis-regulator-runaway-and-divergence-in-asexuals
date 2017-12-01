#include "MoyFitModel6.h"
#include "SelectionModel6.h"
#include <string>
#include <vector>

using namespace std;

double MoyFitModel6(double Ei, vector<double> PopDChrom1, vector<double> PopEChrom1, vector<double> PopSChrom1, vector<double> PopDChrom2, vector<double> PopEChrom2, vector<double> PopSChrom2, double h, int Npop, string FIT, double I, string Modifier)
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
            WMOY += SelectionModel6(Ei,PopDChrom1[indFit],PopEChrom1[indFit],PopSChrom1[indFit],PopDChrom2[indFit],PopEChrom2[indFit],PopSChrom2[indFit],h,I,Modifier);
        }
    }

    WMOY /= Npop;

    return WMOY;
}
