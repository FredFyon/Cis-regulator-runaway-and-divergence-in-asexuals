#include "RecombinaisonAutomixis.h"
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
#include <string>


using namespace std;

Vecteur RecombinaisonAutomixis(double EChrom1, double SChrom1, double EChrom2, double SChrom2, double Rjk, double Rcentro, string AUTOMIXIS)
{
    //ofstream Test("/home/frederic/Bureau/GEFSimulator_CB/bin/Debug/TestR1.txt",ios::app);

    typedef boost::uniform_real<double> UniformDistributionReal;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionReal> UniformGeneratorReal;

    struct timeval Time;
    gettimeofday(&Time,0);
    static RandomGenerator rng(Time.tv_sec+Time.tv_usec*1000);

    UniformDistributionReal LoiUnifReal(0.,1.);
    UniformGeneratorReal TirageUnifReal(rng,LoiUnifReal);

    double EChrom_Gam1,EChrom_Gam2,EChrom_Gam3,EChrom_Gam4,SChrom_Gam1,SChrom_Gam2,SChrom_Gam3,SChrom_Gam4;
    double EChrom1_Desc, SChrom1_Desc, EChrom2_Desc, SChrom2_Desc;

    double ProbaRecombinaison = TirageUnifReal();
    if(ProbaRecombinaison <= Rcentro)
    {
        int Recomb1 = rand() %4 + 1;
        if(Recomb1 == 1)
        {
            EChrom_Gam1 = EChrom2;
            EChrom_Gam2 = EChrom1;
            EChrom_Gam3 = EChrom1;
            EChrom_Gam4 = EChrom2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= Rjk)
            {
                int Recomb2 = rand() %4 + 1;
                if(Recomb2 == 1)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom2;
                }
                else if(Recomb2 == 2)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom2;
                }
                else if(Recomb2 == 3)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom2;
                }
                else if(Recomb2 == 4)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom1;
                }
            }
            else
            {
                SChrom_Gam1 = SChrom2;
                SChrom_Gam2 = SChrom1;
                SChrom_Gam3 = SChrom1;
                SChrom_Gam4 = SChrom2;
            }
        }
        else if(Recomb1 == 2)
        {
            EChrom_Gam1 = EChrom2;
            EChrom_Gam2 = EChrom1;
            EChrom_Gam3 = EChrom2;
            EChrom_Gam4 = EChrom1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= Rjk)
            {
                int Recomb2 = rand() %4 + 1;
                if(Recomb2 == 1)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom1;
                }
                else if(Recomb2 == 2)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom2;
                }
                else if(Recomb2 == 3)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom1;
                }
                else if(Recomb2 == 4)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom1;
                }
            }
            else
            {
                SChrom_Gam1 = SChrom2;
                SChrom_Gam2 = SChrom1;
                SChrom_Gam3 = SChrom2;
                SChrom_Gam4 = SChrom1;
            }
        }
        else if(Recomb1 == 3)
        {
            EChrom_Gam1 = EChrom1;
            EChrom_Gam2 = EChrom2;
            EChrom_Gam3 = EChrom1;
            EChrom_Gam4 = EChrom2;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= Rjk)
            {
                int Recomb2 = rand() %4 + 1;
                if(Recomb2 == 1)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom2;
                }
                else if(Recomb2 == 2)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom1;
                }
                else if(Recomb2 == 3)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom2;
                }
                else if(Recomb2 == 4)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom2;
                }
            }
            else
            {
                SChrom_Gam1 = SChrom1;
                SChrom_Gam2 = SChrom2;
                SChrom_Gam3 = SChrom1;
                SChrom_Gam4 = SChrom2;
            }
        }
        else if(Recomb1 == 4)
        {
            EChrom_Gam1 = EChrom1;
            EChrom_Gam2 = EChrom2;
            EChrom_Gam3 = EChrom2;
            EChrom_Gam4 = EChrom1;

            double ProbaRecombinaison2 = TirageUnifReal();
            if(ProbaRecombinaison2 <= Rjk)
            {
                int Recomb2 = rand() %4 + 1;
                if(Recomb2 == 1)
                {
                    SChrom_Gam1 = SChrom2;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom1;
                    SChrom_Gam4 = SChrom1;
                }
                else if(Recomb2 == 2)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom1;
                }
                else if(Recomb2 == 3)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom2;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom1;
                }
                else if(Recomb2 == 4)
                {
                    SChrom_Gam1 = SChrom1;
                    SChrom_Gam2 = SChrom1;
                    SChrom_Gam3 = SChrom2;
                    SChrom_Gam4 = SChrom2;
                }
            }
            else
            {
                SChrom_Gam1 = SChrom1;
                SChrom_Gam2 = SChrom2;
                SChrom_Gam3 = SChrom2;
                SChrom_Gam4 = SChrom1;
            }
        }
    }
    else
    {
        EChrom_Gam1 = EChrom1;
        EChrom_Gam2 = EChrom1;
        EChrom_Gam3 = EChrom2;
        EChrom_Gam4 = EChrom2;

        double ProbaRecombinaison2 = TirageUnifReal();
        if(ProbaRecombinaison2 <= Rjk)
        {
            int Recomb2 = rand() %4 + 1;
            if(Recomb2 == 1)
            {
                SChrom_Gam1 = SChrom2;
                SChrom_Gam2 = SChrom1;
                SChrom_Gam3 = SChrom1;
                SChrom_Gam4 = SChrom2;
            }
            else if(Recomb2 == 2)
            {
                SChrom_Gam1 = SChrom2;
                SChrom_Gam2 = SChrom1;
                SChrom_Gam3 = SChrom2;
                SChrom_Gam4 = SChrom1;
            }
            else if(Recomb2 == 3)
            {
                SChrom_Gam1 = SChrom1;
                SChrom_Gam2 = SChrom2;
                SChrom_Gam3 = SChrom1;
                SChrom_Gam4 = SChrom2;
            }
            else if(Recomb2 == 4)
            {
                SChrom_Gam1 = SChrom1;
                SChrom_Gam2 = SChrom2;
                SChrom_Gam3 = SChrom2;
                SChrom_Gam4 = SChrom1;
            }
        }
        else
        {
            SChrom_Gam1 = SChrom1;
            SChrom_Gam2 = SChrom1;
            SChrom_Gam3 = SChrom2;
            SChrom_Gam4 = SChrom2;
        }
    }

    if(AUTOMIXIS == "CF")
    {
        int Par = rand() %4 + 1;
        if(Par == 1)
        {
            EChrom1_Desc = EChrom_Gam1;
            SChrom1_Desc = SChrom_Gam1;
            int Par2 = rand() %2 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam3;
                SChrom2_Desc = SChrom_Gam3;
            }
            else
            {
                EChrom2_Desc = EChrom_Gam4;
                SChrom2_Desc = SChrom_Gam4;
            }
        }
        else if(Par == 2)
        {
            EChrom1_Desc = EChrom_Gam2;
            SChrom1_Desc = SChrom_Gam2;
            int Par2 = rand() %2 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam3;
                SChrom2_Desc = SChrom_Gam3;
            }
            else
            {
                EChrom2_Desc = EChrom_Gam4;
                SChrom2_Desc = SChrom_Gam4;
            }
        }
        else if(Par == 3)
        {
            EChrom1_Desc = EChrom_Gam3;
            SChrom1_Desc = SChrom_Gam3;
            int Par2 = rand() %2 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam1;
                SChrom2_Desc = SChrom_Gam1;
            }
            else
            {
                EChrom2_Desc = EChrom_Gam2;
                SChrom2_Desc = SChrom_Gam2;
            }
        }
        else if(Par == 4)
        {
            EChrom1_Desc = EChrom_Gam4;
            SChrom1_Desc = SChrom_Gam4;
            int Par2 = rand() %2 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam1;
                SChrom2_Desc = SChrom_Gam1;
            }
            else
            {
                EChrom2_Desc = EChrom_Gam2;
                SChrom2_Desc = SChrom_Gam2;
            }
        }
    }
    else if(AUTOMIXIS == "TF")
    {
        int Par = rand() %4 + 1;
        if(Par == 1)
        {
            EChrom1_Desc = EChrom_Gam1;
            SChrom1_Desc = SChrom_Gam1;
            EChrom2_Desc = EChrom_Gam2;
            SChrom2_Desc = SChrom_Gam2;
        }
        else if(Par == 2)
        {
            EChrom1_Desc = EChrom_Gam2;
            SChrom1_Desc = SChrom_Gam2;
            EChrom2_Desc = EChrom_Gam1;
            SChrom2_Desc = SChrom_Gam1;
        }
        else if(Par == 3)
        {
            EChrom1_Desc = EChrom_Gam3;
            SChrom1_Desc = SChrom_Gam3;
            EChrom2_Desc = EChrom_Gam4;
            SChrom2_Desc = SChrom_Gam4;
        }
        else if(Par == 4)
        {
            EChrom1_Desc = EChrom_Gam4;
            SChrom1_Desc = SChrom_Gam4;
            EChrom2_Desc = EChrom_Gam3;
            SChrom2_Desc = SChrom_Gam3;
        }
    }
    else if(AUTOMIXIS == "RANDOM")
    {
        int Par = rand() %4 + 1;
        if(Par == 1)
        {
            EChrom1_Desc = EChrom_Gam1;
            SChrom1_Desc = SChrom_Gam1;
            int Par2 = rand() %3 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam2;
                SChrom2_Desc = SChrom_Gam2;
            }
            else if(Par2 == 2)
            {
                EChrom2_Desc = EChrom_Gam3;
                SChrom2_Desc = SChrom_Gam3;
            }
            else if(Par2 == 3)
            {
                EChrom2_Desc = EChrom_Gam4;
                SChrom2_Desc = SChrom_Gam4;
            }
        }
        else if(Par == 2)
        {
            EChrom1_Desc = EChrom_Gam2;
            SChrom1_Desc = SChrom_Gam2;
            int Par2 = rand() %3 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam1;
                SChrom2_Desc = SChrom_Gam1;
            }
            else if(Par2 == 2)
            {
                EChrom2_Desc = EChrom_Gam3;
                SChrom2_Desc = SChrom_Gam3;
            }
            else if(Par2 == 3)
            {
                EChrom2_Desc = EChrom_Gam4;
                SChrom2_Desc = SChrom_Gam4;
            }
        }
        else if(Par == 3)
        {
            EChrom1_Desc = EChrom_Gam3;
            SChrom1_Desc = SChrom_Gam3;
            int Par2 = rand() %3 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam1;
                SChrom2_Desc = SChrom_Gam1;
            }
            else if(Par2 == 2)
            {
                EChrom2_Desc = EChrom_Gam2;
                SChrom2_Desc = SChrom_Gam2;
            }
            else if(Par2 == 3)
            {
                EChrom2_Desc = EChrom_Gam4;
                SChrom2_Desc = SChrom_Gam4;
            }
        }
        else if(Par == 4)
        {
            EChrom1_Desc = EChrom_Gam4;
            SChrom1_Desc = SChrom_Gam4;
            int Par2 = rand() %3 + 1;
            if(Par2 == 1)
            {
                EChrom2_Desc = EChrom_Gam1;
                SChrom2_Desc = SChrom_Gam1;
            }
            else if(Par2 == 2)
            {
                EChrom2_Desc = EChrom_Gam2;
                SChrom2_Desc = SChrom_Gam2;
            }
            else if(Par2 == 3)
            {
                EChrom2_Desc = EChrom_Gam3;
                SChrom2_Desc = SChrom_Gam3;
            }
        }
    }
    else if(AUTOMIXIS == "PostMD")
    {
        int Par = rand() %4 + 1;
        if(Par == 1)
        {
            EChrom1_Desc = EChrom_Gam1;
            SChrom1_Desc = SChrom_Gam1;
            EChrom2_Desc = EChrom_Gam1;
            SChrom2_Desc = SChrom_Gam1;
        }
        else if(Par == 2)
        {
            EChrom1_Desc = EChrom_Gam2;
            SChrom1_Desc = SChrom_Gam2;
            EChrom2_Desc = EChrom_Gam2;
            SChrom2_Desc = SChrom_Gam2;
        }
        else if(Par == 3)
        {
            EChrom1_Desc = EChrom_Gam3;
            SChrom1_Desc = SChrom_Gam3;
            EChrom2_Desc = EChrom_Gam3;
            SChrom2_Desc = SChrom_Gam3;
        }
        else if(Par == 4)
        {
            EChrom1_Desc = EChrom_Gam4;
            SChrom1_Desc = SChrom_Gam4;
            EChrom2_Desc = EChrom_Gam4;
            SChrom2_Desc = SChrom_Gam4;
        }
    }
    else if(AUTOMIXIS == "PreMD")
    {
        double E1C1_0,E1C2_0,S1C1_0,S1C2_0,E2C1_0,E2C2_0,S2C2_0,S2C1_0,E1C1_1,E1C2_1,S1C1_1,S1C2_1,E2C1_1,E2C2_1,S2C2_1,S2C1_1;
        double E1C1_2,E1C2_2,S1C1_2,S1C2_2,E2C1_2,E2C2_2,S2C2_2,S2C1_2,E1C1_3,E1C2_3,S1C1_3,S1C2_3,E2C1_3,E2C2_3,S2C2_3,S2C1_3;

        double ProbaCentro1 = TirageUnifReal();
        if(ProbaCentro1 <= Rcentro)
        {
            double RandRcentro1 = rand() %4 + 1;
            if(RandRcentro1 == 1)
            {
                E1C1_0 = EChrom2;
                E1C1_1 = EChrom1;
                E1C2_0 = EChrom1;
                E1C2_1 = EChrom1;
                S1C1_0 = SChrom2;
                S1C1_1 = SChrom1;
                S1C2_0 = SChrom1;
                S1C2_1 = SChrom1;
            }
            else if(RandRcentro1 == 2)
            {
                E1C1_0 = EChrom1;
                E1C1_1 = EChrom2;
                E1C2_0 = EChrom1;
                E1C2_1 = EChrom1;
                S1C1_0 = SChrom1;
                S1C1_1 = SChrom2;
                S1C2_0 = SChrom1;
                S1C2_1 = SChrom1;
            }
            else if(RandRcentro1 == 3)
            {
                E1C1_0 = EChrom1;
                E1C1_1 = EChrom1;
                E1C2_0 = EChrom2;
                E1C2_1 = EChrom1;
                S1C1_0 = SChrom1;
                S1C1_1 = SChrom1;
                S1C2_0 = SChrom2;
                S1C2_1 = SChrom1;
            }
            else
            {
                E1C1_0 = EChrom1;
                E1C1_1 = EChrom1;
                E1C2_0 = EChrom1;
                E1C2_1 = EChrom2;
                S1C1_0 = SChrom1;
                S1C1_1 = SChrom1;
                S1C2_0 = SChrom1;
                S1C2_1 = SChrom2;
            }
            double RandRcentro2 = rand() %4 + 1;
            if(RandRcentro2 == 1)
            {
                E2C1_0 = EChrom1;
                E2C1_1 = EChrom2;
                E2C2_0 = EChrom2;
                E2C2_1 = EChrom2;
                S2C1_0 = SChrom1;
                S2C1_1 = SChrom2;
                S2C2_0 = SChrom2;
                S2C2_1 = SChrom2;
            }
            else if(RandRcentro2 == 2)
            {
                E2C1_0 = EChrom2;
                E2C1_1 = EChrom1;
                E2C2_0 = EChrom2;
                E2C2_1 = EChrom2;
                S2C1_0 = SChrom2;
                S2C1_1 = SChrom1;
                S2C2_0 = SChrom2;
                S2C2_1 = SChrom2;
            }
            else if(RandRcentro2 == 3)
            {
                E2C1_0 = EChrom2;
                E2C1_1 = EChrom2;
                E2C2_0 = EChrom1;
                E2C2_1 = EChrom2;
                S2C1_0 = SChrom2;
                S2C1_1 = SChrom2;
                S2C2_0 = SChrom1;
                S2C2_1 = SChrom2;
            }
            else
            {
                E2C1_0 = EChrom2;
                E2C1_1 = EChrom2;
                E2C2_0 = EChrom2;
                E2C2_1 = EChrom1;
                S2C1_0 = SChrom2;
                S2C1_1 = SChrom2;
                S2C2_0 = SChrom2;
                S2C2_1 = SChrom1;
            }
        }
        else
        {
            E1C1_0 = EChrom1;
            E1C1_1 = EChrom1;
            E1C2_0 = EChrom1;
            E1C2_1 = EChrom1;
            S1C1_0 = SChrom1;
            S1C1_1 = SChrom1;
            S1C2_0 = SChrom1;
            S1C2_1 = SChrom1;
            E2C1_0 = EChrom2;
            E2C1_1 = EChrom2;
            E2C2_0 = EChrom2;
            E2C2_1 = EChrom2;
            S2C1_0 = SChrom2;
            S2C1_1 = SChrom2;
            S2C2_0 = SChrom2;
            S2C2_1 = SChrom2;
        }

        double ProbaRjk = TirageUnifReal();
        E1C1_2 = E1C1_0;
        E1C1_3 = E1C1_1;
        E1C2_2 = E1C2_0;
        E1C2_3 = E1C2_1;
        E2C1_2 = E2C1_0;
        E2C1_3 = E2C1_1;
        E2C2_2 = E2C2_0;
        E2C2_3 = E2C2_1;

        if(ProbaRjk <= Rjk)
        {
            double RandRjk1 = rand() %16 +1;
            if(RandRjk1 == 1)
            {
                S1C1_2 = S2C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S1C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 2)
            {
                S1C1_2 = S2C1_1;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S1C1_0;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 3)
            {
                S1C1_2 = S2C2_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S1C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 4)
            {
                S1C1_2 = S2C2_1;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S1C1_0;
            }
            else if(RandRjk1 == 5)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S2C1_0;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S1C1_1;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 6)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S2C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S1C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 7)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S2C2_0;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S1C1_1;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 8)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S2C2_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S1C1_1;
            }
            else if(RandRjk1 == 9)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S2C1_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S1C2_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 10)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S2C1_1;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S1C2_0;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 11)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S2C2_0;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S1C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 12)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S2C2_1;
                S1C2_3 = S1C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S1C2_0;
            }
            else if(RandRjk1 == 13)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S2C1_0;
                S2C1_2 = S1C2_1;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 14)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S2C1_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S1C2_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 15)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S2C2_0;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S1C2_1;
                S2C2_3 = S2C2_1;
            }
            else if(RandRjk1 == 16)
            {
                S1C1_2 = S1C1_0;
                S1C1_3 = S1C1_1;
                S1C2_2 = S1C2_0;
                S1C2_3 = S2C2_1;
                S2C1_2 = S2C1_0;
                S2C1_3 = S2C1_1;
                S2C2_2 = S2C2_0;
                S2C2_3 = S1C2_1;
            }
        }
        else
        {
            S1C1_2 = S1C1_0;
            S1C1_3 = S1C1_1;
            S1C2_2 = S1C2_0;
            S1C2_3 = S1C2_1;
            S2C1_2 = S2C1_0;
            S2C1_3 = S2C1_1;
            S2C2_2 = S2C2_0;
            S2C2_3 = S2C2_1;
        }

        double E1,S1,E2,S2;
        double Repro1 = rand() %4 + 1;
        if(Repro1 == 1)
        {
            E1 = E1C1_2;
            S1 = S1C1_2;
        }
        else if(Repro1 == 2)
        {
            E1 = E1C1_3;
            S1 = S1C1_3;
        }
        else if(Repro1 == 3)
        {
            E1 = E1C2_2;
            S1 = S1C2_2;
        }
        else
        {
            E1 = E1C2_3;
            S1 = S1C2_3;
        }

        double Repro2 = rand() %4 + 1;
        if(Repro2 == 1)
        {
            E2 = E2C1_2;
            S2 = S2C1_2;
        }
        else if(Repro2 == 2)
        {
            E2 = E2C1_3;
            S2 = S2C1_3;
        }
        else if(Repro2 == 3)
        {
            E2 = E2C2_2;
            S2 = S2C2_2;
        }
        else
        {
            E2 = E2C2_3;
            S2 = S2C2_3;
        }

        double Stoch = rand() %2 + 1;
        if(Stoch == 1)
        {
            EChrom1_Desc = E1;
            SChrom1_Desc = S1;
            EChrom2_Desc = E2;
            SChrom2_Desc = S2;
        }
        else
        {
            EChrom1_Desc = E2;
            SChrom1_Desc = S2;
            EChrom2_Desc = E1;
            SChrom2_Desc = S1;
        }
    }

    Vecteur V(EChrom1_Desc,SChrom1_Desc,EChrom2_Desc,SChrom2_Desc,0.,0.,0.,0.,0.,0.);

    return V;
}
