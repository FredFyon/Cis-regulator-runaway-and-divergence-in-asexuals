#include "RecombinaisonModel3.h"
#include "Vecteur.h"

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/math/distributions.hpp>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <fstream>


using namespace std;

Vecteur RecombinaisonModel3(double TChrom1_Par1, double EChrom1_Par1, double SChrom1_Par1, double TChrom2_Par1, double EChrom2_Par1, double SChrom2_Par1, double TChrom1_Par2, double EChrom1_Par2, double SChrom1_Par2, double TChrom2_Par2, double EChrom2_Par2, double SChrom2_Par2, double Rij, double Rjk)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/TestR3.txt",ios::app);
    typedef boost::uniform_real<double> UniformDistributionReal;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionReal> UniformGeneratorReal;

    struct timeval Time;
    gettimeofday(&Time,0);
    static RandomGenerator rng(Time.tv_sec+Time.tv_usec*1000);

    UniformDistributionReal LoiUnifReal(0.,1.);
    UniformGeneratorReal TirageUnifReal(rng,LoiUnifReal);

    double TChrom1_Inter, EChrom1_Inter, SChrom1_Inter, TChrom2_Inter, EChrom2_Inter, SChrom2_Inter;
    double TChrom1_Desc, EChrom1_Desc, SChrom1_Desc, TChrom2_Desc, EChrom2_Desc, SChrom2_Desc;

    int chrom1 = rand() %2 + 1;
    int chrom2 = rand() %2 + 1;

    if(chrom1 == 1)
    {
        TChrom1_Inter = TChrom1_Par1;

        double ProbaRecombinaison1 = TirageUnifReal();
        //Test << ProbaRecombinaison1 << endl;
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom1_Inter = EChrom2_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom1_Inter = SChrom1_Par1;
            }
            else
            {
                SChrom1_Inter = SChrom2_Par1;
            }
        }
        else
        {
            EChrom1_Inter = EChrom1_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom1_Inter = SChrom2_Par1;
            }
            else
            {
                SChrom1_Inter = SChrom1_Par1;
            }
        }
    }
    else
    {
        TChrom1_Inter = TChrom2_Par1;

        double ProbaRecombinaison1 = TirageUnifReal();
        //Test << ProbaRecombinaison1 << endl;
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom1_Inter = EChrom1_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom1_Inter = SChrom2_Par1;
            }
            else
            {
                SChrom1_Inter = SChrom1_Par1;
            }
        }
        else
        {
            EChrom1_Inter = EChrom2_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom1_Inter = SChrom1_Par1;
            }
            else
            {
                SChrom1_Inter = SChrom2_Par1;
            }
        }
    }

    if(chrom2 == 1)
    {
        TChrom2_Inter = TChrom1_Par2;

        double ProbaRecombinaison1 = TirageUnifReal();
        //Test << ProbaRecombinaison1 << endl;
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom2_Inter = EChrom2_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom2_Inter = SChrom1_Par2;
            }
            else
            {
                SChrom2_Inter = SChrom2_Par2;
            }
        }
        else
        {
            EChrom2_Inter = EChrom1_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom2_Inter = SChrom2_Par2;
            }
            else
            {
                SChrom2_Inter = SChrom1_Par2;
            }
        }
    }
    else
    {
        TChrom2_Inter = TChrom2_Par2;

        double ProbaRecombinaison1 = TirageUnifReal();
        //Test << ProbaRecombinaison1 << endl;
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom2_Inter = EChrom1_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom2_Inter = SChrom2_Par2;
            }
            else
            {
                SChrom2_Inter = SChrom1_Par2;
            }
        }
        else
        {
            EChrom2_Inter = EChrom2_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            //Test << ProbaRecombinaison2 << endl;
            if(ProbaRecombinaison2 <= Rjk)
            {
                SChrom2_Inter = SChrom1_Par2;
            }
            else
            {
                SChrom2_Inter = SChrom2_Par2;
            }
        }
    }

    int Par = rand() %2 + 1;
    if(Par == 1)
    {
        TChrom1_Desc = TChrom1_Inter;
        EChrom1_Desc = EChrom1_Inter;
        SChrom1_Desc = SChrom1_Inter;
        TChrom2_Desc = TChrom2_Inter;
        EChrom2_Desc = EChrom2_Inter;
        SChrom2_Desc = SChrom2_Inter;
    }
    else
    {
        TChrom1_Desc = TChrom2_Inter;
        EChrom1_Desc = EChrom2_Inter;
        SChrom1_Desc = SChrom2_Inter;
        TChrom2_Desc = TChrom1_Inter;
        EChrom2_Desc = EChrom1_Inter;
        SChrom2_Desc = SChrom1_Inter;
    }

    Vecteur V(TChrom1_Desc,EChrom1_Desc,SChrom1_Desc,TChrom2_Desc,EChrom2_Desc,SChrom2_Desc,0.,0.,0.,0.);

    return V;
}
