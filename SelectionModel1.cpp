#include "SelectionModel1.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "Vecteur.h"

using namespace std;

Vecteur SelectionModel1(double EChrom1, double SChrom1, double EChrom2, double SChrom2, double h)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/h.txt",ios::app);
    double S = max(SChrom1,SChrom2);

    double y,W;

    double H;

    if((EChrom1 == 0) && (EChrom2 == 0))
    {
        cout << "Promoters attractivity for transcription factors has collapsed. Model can't calculate individual fitness." << endl;
    }
    else
    {
        if(S == SChrom1)
        {
            y = EChrom1/(EChrom1+EChrom2);
            W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
            H = pow(y,-log(h)/log(2));
            //Test << pow(y,-log(h)/log(2)) << endl;
        }
        else
        {
            y = EChrom2/(EChrom1+EChrom2);
            W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
            H = pow(y,-log(h)/log(2));
            //Test << pow(y,-log(h)/log(2)) << endl;
        }
    }

    //Test << SChrom1 << " " << SChrom2 << " " << W << endl;

    //Test << W << endl;

    //cout << EChrom1 << " " << EChrom2 << " " << SChrom1 << " " << SChrom2 << " " << W << endl;

    Vecteur Resultats(W,H,0.,0.,0.,0.,0.,0.,0.,0.);

    return Resultats;
}

