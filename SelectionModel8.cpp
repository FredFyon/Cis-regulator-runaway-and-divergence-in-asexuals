#include "SelectionModel8.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include "Vecteur.h"
#include <stdlib.h>
#include <string>

using namespace std;

Vecteur SelectionModel8(double Ei, double SexChrom1, double TChrom1, double EChrom1, double SChrom1, double SexChrom2, double TChrom2, double EChrom2, double SChrom2, double h, double I, string SEX)
{
    //cout << TChrom1 << " " << EChrom1 << " " << SChrom1 << " " << TChrom2 << " " << EChrom2 << " " << SChrom2 << endl;

    srand (time(NULL));

    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/h.txt",ios::app);
    double S = max(SChrom1,SChrom2);

    double Male(0),Female(0);

    double y,W,E;

    double H;

    double Xmale(0),Y(0),Xfem1(0),Xfem2(0);

    if((EChrom1 == 0) && (EChrom2 == 0))
    {
        cout << "Promoters attractivity for transcription factors has collapsed. Model can't calculate individual fitness." << endl;
    }
    else
    {
        if(SEX == "1")
        {
            if(SexChrom1 != SexChrom2)
            {
                Male = 1;
                if(SexChrom1 == 0)
                {
                    if(S = SChrom1)
                    {
                        y = ((log10(TChrom1)+log10(TChrom2))/2.)*EChrom1/(((log10(TChrom1)+log10(TChrom2))/2)*EChrom1+EChrom2);
                        W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
                    }
                    else
                    {
                        y = EChrom2/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom1+EChrom2);
                        W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
                    }
                    //E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1) + log10(EChrom2);
                    E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1) + log10(EChrom2);
                    Xmale = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1);
                    Y = log10(EChrom2);
                    H = pow(EChrom2/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom1+EChrom2),-log(h)/log(2));
                }
                else
                {
                    if(S = SChrom1)
                    {
                        y = EChrom1/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2+EChrom1);
                        W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
                    }
                    else
                    {
                        y = ((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2+EChrom1);
                        W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
                    }
                    //E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2) + log10(EChrom1);
                    E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2) + log10(EChrom1);
                    Xmale = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2);
                    Y = log10(EChrom1);
                    H = pow(EChrom1/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2+EChrom1),-log(h)/log(2));
                }
            }
            else
            {
                Female = 1;
                if(S = SChrom1)
                {
                    y = EChrom1/(EChrom1+EChrom2);
                    W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
                }
                else
                {
                    y = EChrom2/(EChrom1+EChrom2);
                    W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
                }
                //E = log10(EChrom2) + log10(EChrom1);
                E = log10(EChrom2) + log10(EChrom1);
                Xfem1 = log10(EChrom1);
                Xfem2 = log10(EChrom2);
                H = 0;
            }
        }

        else if(SEX == "2")
        {
            if(S == SChrom1)
            {
                y = EChrom1/(EChrom1+EChrom2);
                W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
            }
            else
            {
                y = EChrom2/(EChrom1+EChrom2);
                W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
            }

            if(SexChrom1 == SexChrom2)
            {
                Female = 1;
                //E = ((log10(TChrom1)+log10(TChrom2))/2.)*(log10(EChrom1) + log10(EChrom2));
                E = ((log10(TChrom1)+log10(TChrom2))/2.)*(log10(EChrom1) + log10(EChrom2));
                Xfem1 = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1);
                Xfem2 = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2);
                H = 0;
            }
            else
            {
                Male = 1;
                //E = log10(EChrom1) + log10(EChrom2);
                E = log10(EChrom1) + log10(EChrom2);
                if(SexChrom1 == 0)
                {
                    Xmale = log10(EChrom1);
                    Y = log10(EChrom2);
                    H = pow(EChrom2/(EChrom1+EChrom2),-log(h)/log(2));
                }
                else
                {
                    Xmale = log10(EChrom2);
                    Y = log10(EChrom1);
                    H = pow(EChrom1/(EChrom1+EChrom2),-log(h)/log(2));
                }
            }
        }

        else if(SEX == "3")
        {
            if(SexChrom1 == SexChrom2)
            {
                Female = 1;
                double choix = rand() %2 + 1;
                if(choix == 1)
                {
                    if(S == SChrom1)
                    {
                        y = ((log10(TChrom1)+log10(TChrom2))/2.)*EChrom1/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom1+EChrom2);
                        W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
                    }
                    else
                    {
                        y = EChrom2/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom1+EChrom2);
                        W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
                    }
                    //E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1) + log10(EChrom2);
                    E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1) + log10(EChrom2);
                    Xfem1 = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom1);
                    Xfem2 = log10(EChrom2);
                }
                else
                {
                    if(S == SChrom1)
                    {
                        y = EChrom1/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2+EChrom1);
                        W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
                    }
                    else
                    {
                        y = ((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2/(((log10(TChrom1)+log10(TChrom2))/2.)*EChrom2+EChrom1);
                        W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
                    }
                    //E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2) + log10(EChrom1);
                    E = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2) + log10(EChrom1);
                    Xfem1 = ((log10(TChrom1)+log10(TChrom2))/2.)*log10(EChrom2);
                    Xfem2 = log10(EChrom1);
                }
                H = 0;
            }
            else
            {
                Male = 1;
                if(S == SChrom1)
                {
                    y = EChrom1/(EChrom1+EChrom2);
                    W = 1 - SChrom2 + (SChrom2 - SChrom1)*pow(y,-log(h)/log(2));
                }
                else
                {
                    y = EChrom2/(EChrom1+EChrom2);
                    W = 1 - SChrom1 + (SChrom1 - SChrom2)*pow(y,-log(h)/log(2));
                }
                //E = log10(EChrom2) + log10(EChrom1);
                E = log10(EChrom2) + log10(EChrom1);
                if(SexChrom1 == 0)
                {
                    Xmale = log10(EChrom1);
                    Y = log10(EChrom2);
                    H = pow(EChrom2/(EChrom1+EChrom2),-log(h)/log(2));
                }
                else
                {
                    Xmale = log10(EChrom2);
                    Y = log10(EChrom1);
                    H = pow(EChrom1/(EChrom1+EChrom2),-log(h)/log(2));
                }
            }
        }
    }

    //cout << Xmale << " " << Y << " " << Xfem1 << " " << Xfem2 << endl;

    double Wconstraint;

    Wconstraint = exp(-I*pow(E-2*log10(Ei),2));

    double Wtot = W*Wconstraint;

    //cout << E << endl;

    //Test << W << endl;

    //cout << EChrom1 << " " << EChrom2 << " " << SChrom1 << " " << SChrom2 << " " << W << endl;

    Vecteur Resultats(Wtot,H,Xmale,Y,Xfem1,Xfem2,Male,Female,0.,0.);

    return Resultats;
}
