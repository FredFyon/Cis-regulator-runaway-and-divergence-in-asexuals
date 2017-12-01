#include "SelectionModel5.h"
#include <string>
#include <cmath>
#include <iostream>
#include "math.h"
#include "Vecteur.h"

using namespace std;

Vecteur SelectionModel5(double Ei, double DChrom1, double EChrom1, double SChrom1, double DChrom2, double EChrom2, double SChrom2, double h, double I)
{
    double S = max(SChrom1,SChrom2);
    /*if(DChrom1 < 0)
    {
        DChrom1 = 0;
    }
    else if(DChrom1 > 1)
    {
        DChrom1 = 1;
    }

    if(DChrom2 < 0)
    {
        DChrom2 = 0;
    }
    else if(DChrom2 > 1)
    {
        DChrom2 = 1;
    }*/

    double y,Wes,H,hh;
    double alpha = (log10(DChrom1)+log10(DChrom2))/2.;

    if((EChrom1 == 0) && (EChrom2 == 0))
    {
        cout << "Promoters attractivity for transcription factors has collapsed. Model can't calculate individual fitness." << endl;
    }
    else
    {
        if(S == SChrom1)
        {
            y = EChrom1/(EChrom1+EChrom2);
            if(SChrom1 != SChrom2)
            {
                if(EChrom1 == EChrom2)
                {
                    hh = h;
                }
                else
                {
                    if(alpha > 1)
                    {
                        hh = h;
                    }
                    else
                    {
                        hh = (h + (pow(y,-log(h)/log(2))-h)*(1-alpha));
                    }
                }
            }
            else
            {
                hh = 0.;
            }

            Wes = 1 - SChrom2 + (SChrom2 - SChrom1)*hh;
            //cout << Wes << endl;
        }
        else
        {
            y = EChrom2/(EChrom1+EChrom2);
            if(SChrom1 != SChrom2)
            {
                if(EChrom1 == EChrom2)
                {
                    hh = h;
                }
                else
                {
                    if(alpha > 1)
                    {
                        hh = h;
                    }
                    else
                    {
                        hh = (h + (pow(y,-log(h)/log(2))-h)*(1-alpha));
                    }
                }
            }
            else
            {
                hh = 0.;
            }

            Wes = 1 - SChrom1 + (SChrom1 - SChrom2)*hh;
        }
    }

    double E = log10(EChrom1) + log10(EChrom2);
    double Wconstraint;

    Wconstraint = exp(-I*pow(E-2*log10(Ei),2));

    //cout << Wconstraint << endl;

    double W = Wes*Wconstraint;

    Vecteur Res(W,hh,0.,0.,0.,0.,0.,0.,0.,0.);

    //cout << hh << " " << SChrom1 << " " << SChrom2 << endl;

    //cout << W << endl;

    return Res;
}
