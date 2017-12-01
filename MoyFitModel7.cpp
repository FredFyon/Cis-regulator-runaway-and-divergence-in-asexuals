#include "MoyFitModel7.h"
#include "SelectionModel7.h"
#include <string>
#include <vector>

using namespace std;

double MoyFitModel7(double Ei, vector<double> PopE1Chrom1, vector<double> PopE2Chrom1, vector<double> PopSChrom1, vector<double> PopE1Chrom2, vector<double> PopE2Chrom2, vector<double> PopSChrom2, double h, int Npop, string FIT, double I)
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
            WMOY += SelectionModel7(Ei,PopE1Chrom1[indFit],PopE2Chrom1[indFit],PopSChrom1[indFit],PopE1Chrom2[indFit],PopE2Chrom2[indFit],PopSChrom2[indFit],h,I);
        }
    }

    WMOY /= Npop;

    return WMOY;
}
