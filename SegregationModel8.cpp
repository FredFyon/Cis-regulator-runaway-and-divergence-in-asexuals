#include "SegregationModel8.h"
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

Vecteur SegregationModel8(double SexChrom1_Par1,double TChrom1_Par1,double EChrom1_Par1, double SChrom1_Par1, double SexChrom2_Par1,double TChrom2_Par1,double EChrom2_Par1, double SChrom2_Par1, double SexChrom1_Par2, double TChrom1_Par2,double EChrom1_Par2, double SChrom1_Par2, double SexChrom2_Par2, double TChrom2_Par2,double EChrom2_Par2, double SChrom2_Par2, double Rte, double R)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/TestR1.txt",ios::app);

    //cout << TChrom1_Par2 << " " << EChrom1_Par2 << " " << SChrom1_Par2 << " " << TChrom2_Par2 << " " << EChrom2_Par2 << " " << SChrom2_Par2 << endl;

    typedef boost::uniform_real<double> UniformDistributionReal;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionReal> UniformGeneratorReal;

    struct timeval Time;
    gettimeofday(&Time,0);
    static RandomGenerator rng(Time.tv_sec+Time.tv_usec*1000);

    UniformDistributionReal LoiUnifReal(0.,1.);
    UniformGeneratorReal TirUnifReal(rng,LoiUnifReal);

    double SexChrom1_Inter, TChrom1_Inter, EChrom1_Inter, SChrom1_Inter, SexChrom2_Inter, TChrom2_Inter, EChrom2_Inter, SChrom2_Inter,SexChrom1_Desc, TChrom1_Desc, EChrom1_Desc, SChrom1_Desc, SexChrom2_Desc, TChrom2_Desc, EChrom2_Desc, SChrom2_Desc;

    int Par1 = rand() %2 + 1;
    if(Par1 == 1)
    {
        SexChrom1_Inter = SexChrom1_Par1;
        TChrom1_Inter = TChrom1_Par1;
        if(SexChrom1_Par1 == SexChrom2_Par1)
        {
            double ProbaRecomb = TirUnifReal();
            if(ProbaRecomb < Rte)
            {
                EChrom1_Inter = EChrom2_Par1;
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
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
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
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
            EChrom1_Inter = EChrom1_Par1;
            SChrom1_Inter = SChrom1_Par1;
        }
    }
    else
    {
        SexChrom1_Inter = SexChrom2_Par1;
        TChrom1_Inter = TChrom2_Par1;
        if(SexChrom1_Par1 == SexChrom2_Par1)
        {
            double ProbaRecomb = TirUnifReal();
            if(ProbaRecomb < Rte)
            {
                EChrom1_Inter = EChrom1_Par1;
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
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
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
                {
                    SChrom1_Inter = SChrom1_Par1;
                }
                else
                {
                    SChrom1_Inter = SChrom2_Par1;
                }
            }
        }
        else
        {
            EChrom1_Inter = EChrom2_Par1;
            SChrom1_Inter = SChrom2_Par1;
        }
    }

    int Par2 = rand() %2 + 1;
    if(Par2 == 1)
    {
        SexChrom2_Inter = SexChrom1_Par2;
        TChrom2_Inter = TChrom1_Par2;
        if(SexChrom1_Par2 == SexChrom2_Par2)
        {
            double ProbaRecomb = TirUnifReal();
            if(ProbaRecomb < Rte)
            {
                EChrom2_Inter = EChrom2_Par2;
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
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
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
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
            EChrom2_Inter = EChrom1_Par2;
            SChrom2_Inter = SChrom1_Par2;
        }
    }
    else
    {
        SexChrom2_Inter = SexChrom2_Par2;
        TChrom2_Inter = TChrom2_Par2;
        if(SexChrom1_Par2 == SexChrom2_Par2)
        {
            double ProbaRecomb = TirUnifReal();
            if(ProbaRecomb < Rte)
            {
                EChrom2_Inter = EChrom1_Par2;
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
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
                double ProbaRecomb2 = TirUnifReal();
                if(ProbaRecomb2 < R)
                {
                    SChrom2_Inter = SChrom1_Par2;
                }
                else
                {
                    SChrom2_Inter = SChrom2_Par2;
                }
            }
        }
        else
        {
            EChrom2_Inter = EChrom2_Par2;
            SChrom2_Inter = SChrom2_Par2;
        }
    }

    int Par = rand() %2 + 1;
    if(Par == 1)
    {
        SexChrom1_Desc = SexChrom1_Inter;
        TChrom1_Desc = TChrom1_Inter;
        EChrom1_Desc = EChrom1_Inter;
        SChrom1_Desc = SChrom1_Inter;
        SexChrom2_Desc = SexChrom2_Inter;
        TChrom2_Desc = TChrom2_Inter;
        EChrom2_Desc = EChrom2_Inter;
        SChrom2_Desc = SChrom2_Inter;
    }
    else
    {
        SexChrom1_Desc = SexChrom2_Inter;
        TChrom1_Desc = TChrom2_Inter;
        EChrom1_Desc = EChrom2_Inter;
        SChrom1_Desc = SChrom2_Inter;
        SexChrom2_Desc = SexChrom1_Inter;
        TChrom2_Desc = TChrom1_Inter;
        EChrom2_Desc = EChrom1_Inter;
        SChrom2_Desc = SChrom1_Inter;
    }

    //cout << TChrom1_Inter << " " << EChrom1_Inter << " " << SChrom1_Inter << " " << TChrom2_Inter << " " << EChrom2_Inter << " " << SChrom2_Inter << endl;
    //cout << TChrom1_Desc << " " << EChrom1_Desc << " " << SChrom1_Desc << " " << TChrom2_Desc << " " << EChrom2_Desc << " " << SChrom2_Desc << endl;

    Vecteur V(SexChrom1_Desc,TChrom1_Desc,EChrom1_Desc,SChrom1_Desc,SexChrom2_Desc,TChrom2_Desc,EChrom2_Desc,SChrom2_Desc,0.,0.);

    return V;
}
