#include "MoyFitModel3.h"
#include "SelectionModel3.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

double MoyFitModel3(vector<double> PopTChrom1, vector<double> PopEChrom1, vector<double> PopSChrom1, vector<double> PopTChrom2, vector<double> PopEChrom2, vector<double> PopSChrom2, double h, int Npop, string FIT, double I, double Concaveness, string SHAPE)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/TestW3.txt",ios::app);

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
            double a = SelectionModel3(PopTChrom1[indFit],PopEChrom1[indFit],PopSChrom1[indFit],PopTChrom2[indFit],PopEChrom2[indFit],PopSChrom2[indFit],h,I,Concaveness,SHAPE);
            WMOY += a;
            //cout << indFit << " " << a << endl;
        }
    }

    double WMoy = WMOY/Npop;

    //Test << WMoy << endl;

    return WMoy;
}
