#include "MoyFitModel1.h"
#include "SelectionModel1.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "Vecteur.h"

using namespace std;

double MoyFitModel1(vector<double> EChrom1, vector<double> SChrom1, vector<double> EChrom2, vector<double> SChrom2, double h, double Npop, string FIT)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/TestW1.txt",ios::app);
    double WMOY(0);
    int indFit(0);

    for(indFit = 0 ; indFit < Npop ; ++indFit)
    {
        if(FIT == "Derive")
        {
            WMOY += 1;
        }
        else if(FIT == "Selection")
        {
            WMOY += SelectionModel1(EChrom1[indFit],SChrom1[indFit],EChrom2[indFit],SChrom2[indFit],h).getx1();
        }
    }

    double WMoy = WMOY/Npop;

   // Test << WMoy << endl;


    //cout << WMOY << endl;

    return WMoy;
}
