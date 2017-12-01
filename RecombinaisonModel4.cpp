#include "RecombinaisonModel4.h"
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


using namespace std;

Vecteur RecombinaisonModel4(double RChrom1_Par1, double EChrom1_Par1, double SChrom1_Par1, double RChrom2_Par1, double EChrom2_Par1, double SChrom2_Par1, double RChrom1_Par2, double EChrom1_Par2, double SChrom1_Par2, double RChrom2_Par2, double EChrom2_Par2, double SChrom2_Par2, double Rij)
{
    typedef boost::uniform_real<double> UniformDistributionReal;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionReal> UniformGeneratorReal;

    struct timeval Time;
    gettimeofday(&Time,0);
    static RandomGenerator rng(Time.tv_sec+Time.tv_usec*1000);

    UniformDistributionReal LoiUnifReal(0.,1.);
    UniformGeneratorReal TirageUnifReal(rng,LoiUnifReal);

    double RChrom1_Inter, EChrom1_Inter, SChrom1_Inter, RChrom2_Inter, EChrom2_Inter, SChrom2_Inter;
    double RChrom1_Desc, EChrom1_Desc, SChrom1_Desc, RChrom2_Desc, EChrom2_Desc, SChrom2_Desc;

    int chrom1 = rand() %2 + 1;
    int chrom2 = rand() %2 + 1;

    double ProbaConversion1 = TirageUnifReal();
    if(ProbaConversion1 <= (RChrom1_Par1+RChrom2_Par1)/2.)
    {
        int conv = rand() %2 + 1;
        if(conv == 1)
        {
            EChrom2_Par1 = EChrom1_Par1;
        }
        else
        {
            EChrom1_Par1 = EChrom2_Par1;
        }
    }

    double ProbaConversion2 = TirageUnifReal();
    if(ProbaConversion2 <= (RChrom1_Par2+RChrom2_Par2)/2.)
    {
        int conv = rand() %2 + 1;
        if(conv == 1)
        {
            EChrom2_Par2 = EChrom1_Par2;
        }
        else
        {
            EChrom1_Par2 = EChrom2_Par2;
        }
    }

    if(chrom1 == 1)
    {
        RChrom1_Inter = RChrom1_Par1;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom1_Inter = EChrom2_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= (RChrom1_Par1+RChrom2_Par1)/2.)
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
            if(ProbaRecombinaison2 <= (RChrom1_Par1+RChrom2_Par1)/2.)
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
        RChrom1_Inter = RChrom2_Par1;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom1_Inter = EChrom1_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= (RChrom1_Par1+RChrom2_Par1)/2.)
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
            if(ProbaRecombinaison2 <= (RChrom1_Par1+RChrom2_Par1)/2.)
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
        RChrom2_Inter = RChrom1_Par2;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom2_Inter = EChrom2_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= (RChrom1_Par2+RChrom2_Par2)/2.)
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
            if(ProbaRecombinaison2 <= (RChrom1_Par2+RChrom2_Par2)/2.)
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
        RChrom2_Inter = RChrom2_Par2;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rij)
        {
            EChrom2_Inter = EChrom1_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= (RChrom1_Par2+RChrom2_Par2)/2.)
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
            if(ProbaRecombinaison2 <= (RChrom1_Par2+RChrom2_Par2)/2.)
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
        RChrom1_Desc = RChrom1_Inter;
        EChrom1_Desc = EChrom1_Inter;
        SChrom1_Desc = SChrom1_Inter;
        RChrom2_Desc = RChrom2_Inter;
        EChrom2_Desc = EChrom2_Inter;
        SChrom2_Desc = SChrom2_Inter;
    }
    else
    {
        RChrom1_Desc = RChrom2_Inter;
        EChrom1_Desc = EChrom2_Inter;
        SChrom1_Desc = SChrom2_Inter;
        RChrom2_Desc = RChrom1_Inter;
        EChrom2_Desc = EChrom1_Inter;
        SChrom2_Desc = SChrom1_Inter;
    }

    Vecteur V(RChrom1_Desc,EChrom1_Desc,SChrom1_Desc,RChrom2_Desc,EChrom2_Desc,SChrom2_Desc,0.,0.,0.,0.);

    return V;
}
