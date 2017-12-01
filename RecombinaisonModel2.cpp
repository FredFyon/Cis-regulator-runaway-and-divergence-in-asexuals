#include "RecombinaisonModel2.h"
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

Vecteur RecombinaisonModel2(double EChrom1a_Par1, double SChrom1a_Par1, double EChrom2a_Par1, double SChrom2a_Par1, double EChrom1b_Par1, double SChrom1b_Par1, double EChrom2b_Par1, double SChrom2b_Par1, double EChrom1a_Par2, double SChrom1a_Par2, double EChrom2a_Par2, double SChrom2a_Par2, double EChrom1b_Par2, double SChrom1b_Par2, double EChrom2b_Par2, double SChrom2b_Par2, double Rjk)
{
    typedef boost::uniform_real<double> UniformDistributionReal;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionReal> UniformGeneratorReal;

    struct timeval Time;
    gettimeofday(&Time,0);
    static RandomGenerator rng(Time.tv_sec+Time.tv_usec*1000);

    UniformDistributionReal LoiUnifReal(0.,1.);
    UniformGeneratorReal TirageUnifReal(rng,LoiUnifReal);

    double EChrom1a_Inter,SChrom1a_Inter,EChrom2a_Inter,SChrom2a_Inter,EChrom1b_Inter,SChrom1b_Inter,EChrom2b_Inter,SChrom2b_Inter;
    double EChrom1a_Desc,SChrom1a_Desc,EChrom2a_Desc,SChrom2a_Desc,EChrom1b_Desc,SChrom1b_Desc,EChrom2b_Desc,SChrom2b_Desc;

    int chrom1 = rand() %2 + 1;
    int chrom2 = rand() %2 + 1;

    if(chrom1 == 1)
    {
        EChrom1a_Inter = EChrom1a_Par1;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rjk)
        {
            SChrom1a_Inter = SChrom2a_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom1b_Inter = EChrom1b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
            }
            else
            {
                EChrom1b_Inter = EChrom2b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
            }
        }
        else
        {
            SChrom1a_Inter = SChrom1a_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom1b_Inter = EChrom2b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
            }
            else
            {
                EChrom1b_Inter = EChrom1b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
            }
        }
    }
    else
    {
        EChrom1a_Inter = EChrom2a_Par1;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rjk)
        {
            SChrom1a_Inter = SChrom1a_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom1b_Inter = EChrom2b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
            }
            else
            {
                EChrom1b_Inter = EChrom1b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
            }
        }
        else
        {
            SChrom1a_Inter = SChrom2a_Par1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom1b_Inter = EChrom1b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
            }
            else
            {
                EChrom1b_Inter = EChrom2b_Par1;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom1b_Inter = SChrom1b_Par1;
                }
                else
                {
                    SChrom1b_Inter = SChrom2b_Par1;
                }
            }
        }
    }

    if(chrom2 == 1)
    {
        EChrom2a_Inter = EChrom1a_Par2;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rjk)
        {
            SChrom2a_Inter = SChrom2a_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom2b_Inter = EChrom1b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
            }
            else
            {
                EChrom2b_Inter = EChrom2b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
            }
        }
        else
        {
            SChrom2a_Inter = SChrom1a_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom2b_Inter = EChrom2b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
            }
            else
            {
                EChrom2b_Inter = EChrom1b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
            }
        }
    }
    else
    {
        EChrom2a_Inter = EChrom2a_Par2;

        double ProbaRecombinaison1 = TirageUnifReal();
        if(ProbaRecombinaison1 <= Rjk)
        {
            SChrom2a_Inter = SChrom1a_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom2b_Inter = EChrom2b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom2b_Par2;
;                }
            }
            else
            {
                EChrom2b_Inter = EChrom1b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
            }
        }
        else
        {
            SChrom2a_Inter = SChrom2a_Par2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= 0.5)
            {
                EChrom2b_Inter = EChrom1b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
            }
            else
            {
                EChrom2b_Inter = EChrom2b_Par2;

                double ProbaRecombinaison3 = TirageUnifReal();
                if(ProbaRecombinaison3 <= Rjk)
                {
                    SChrom2b_Inter = SChrom1b_Par2;
                }
                else
                {
                    SChrom2b_Inter = SChrom2b_Par2;
                }
            }
        }
    }

    int Par = rand() %2 + 1;
    if(Par == 1)
    {
        EChrom1a_Desc = EChrom1a_Inter;
        SChrom1a_Desc = SChrom1a_Inter;
        EChrom2a_Desc = EChrom2a_Inter;
        SChrom2a_Desc = SChrom2a_Inter;
        EChrom1b_Desc = EChrom1b_Inter;
        SChrom1b_Desc = SChrom1b_Inter;
        EChrom2b_Desc = EChrom2b_Inter;
        SChrom2b_Desc = SChrom2b_Inter;
    }
    else
    {
        EChrom1a_Desc = EChrom2a_Inter;
        SChrom1a_Desc = SChrom2a_Inter;
        EChrom2a_Desc = EChrom1a_Inter;
        SChrom2a_Desc = SChrom1a_Inter;
        EChrom1b_Desc = EChrom2b_Inter;
        SChrom1b_Desc = SChrom2b_Inter;
        EChrom2b_Desc = EChrom1b_Inter;
        SChrom2b_Desc = SChrom1b_Inter;
    }

    Vecteur V(EChrom1a_Desc,SChrom1a_Desc,EChrom2a_Desc,SChrom2a_Desc,EChrom1b_Desc,SChrom1b_Desc,EChrom2b_Desc,SChrom2b_Desc,0.,0.);

    return V;
}
