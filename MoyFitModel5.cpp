#include "MoyFitModel5.h"
#include "SelectionModel5.h"
#include <string>
#include <vector>
#include "Vecteur.h"

using namespace std;

double MoyFitModel5(double Ei, vector<double> PopDChrom1, vector<double> PopEChrom1, vector<double> PopSChrom1, vector<double> PopDChrom2, vector<double> PopEChrom2, vector<double> PopSChrom2, double h, int Npop, string FIT, double I)
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
            WMOY += SelectionModel5(Ei,PopDChrom1[indFit],PopEChrom1[indFit],PopSChrom1[indFit],PopDChrom2[indFit],PopEChrom2[indFit],PopSChrom2[indFit],h,I).getx1();
        }
    }

    WMOY /= Npop;

    return WMOY;
}
