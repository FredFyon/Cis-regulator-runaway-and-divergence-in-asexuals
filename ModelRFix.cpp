#include "ModelRFix.h"
#include "Matrice.h"
#include "CycleModel4.h"
#include <string>
#include <cstring>

#include <fstream>
#include <iostream>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>


void ModelRFix(double Ei, double mutE, double mutA, double SigR, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, int Nit, int TpsDesinit, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway)
{
    size_t TailleLineEdit = OutputPathway.size() + 1;
    char* FilePathway = new char[TailleLineEdit];
    strncpy(FilePathway,OutputPathway.c_str(),TailleLineEdit);                                  //Construction du fichier de sortie
    ofstream Results(FilePathway,ios::app);

    typedef boost::uniform_int<int> UniformDistributionInt;
    typedef boost::mt19937 RandomGenerator;
    typedef boost::variate_generator<RandomGenerator&, UniformDistributionInt> UniformGeneratorInt;
    static RandomGenerator rng(clock());
    UniformDistributionInt LoiUnifInt(0,Npop-1);
    UniformGeneratorInt TirUnifInt(rng,LoiUnifInt);

    vector<double> RChrom1(Npop),SChrom1(Npop),EChrom1(Npop),RChrom2(Npop),SChrom2(Npop),EChrom2(Npop);
    vector<double> Empty;
    int InitPop(0);
    for(InitPop = 0 ; InitPop < Npop ; ++InitPop)
    {
        RChrom1[InitPop] = Rjk;
        EChrom1[InitPop] = Ei;
        SChrom1[InitPop] = 0.;
        RChrom2[InitPop] = Rjk;
        EChrom2[InitPop] = Ei;
        SChrom2[InitPop] = 0.;
    }

    Matrice PopInit(RChrom1,EChrom1,SChrom1,RChrom2,EChrom2,SChrom2,Empty,Empty,Empty,Empty);

    int ComptFix(0);
    int IndFix(1);
    for(IndFix = 1 ; IndFix <= Nit ; ++IndFix)
    {
        if(IndFix %(Nit/10) == 0)
        {
            Results << "Fixation " << 10 - (Nit-IndFix)/(Nit/10) << " en cours" << endl;
        }
        Matrice Pop(Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty,Empty);
        int IndDesinit(0);
        for(IndDesinit = 0 ; IndDesinit < TpsDesinit ; ++IndDesinit)
        {
            if(IndDesinit == 0)
            {
                Pop = PopInit;
            }
            Pop = CycleModel4(Pop,Ei,0.,mutE,mutA,SigR,SigE,h,S,Rij,Self,I,Npop,FIT,ALLELES);
        }

        //cout << "desinit done" << endl;

        vector<double> RChrom1Fix = Pop.getVec1();
        vector<double> EChrom1Fix = Pop.getVec2();
        vector<double> SChrom1Fix = Pop.getVec3();
        vector<double> RChrom2Fix = Pop.getVec4();
        vector<double> EChrom2Fix = Pop.getVec5();
        vector<double> SChrom2Fix = Pop.getVec6();

        double chrom = rand() %2 + 1;
        double Mutant = TirUnifInt();
        if(chrom == 1)
        {
            RChrom1Fix[Mutant] = RChrom1[Mutant] + SigR;
        }
        else if(chrom == 2)
        {
            RChrom2Fix[Mutant] = RChrom2[Mutant] + SigR;
        }

        //cout << "Mutation done on " << Mutant << endl;

        Matrice PopFix(RChrom1Fix,EChrom1Fix,SChrom1Fix,RChrom2Fix,EChrom2Fix,SChrom2Fix,Empty,Empty,Empty,Empty);

        int NbAllMut;

        do
        {
            NbAllMut = 0;
            PopFix = CycleModel4(PopFix,Ei,0.,mutE,mutA,SigR,SigE,h,S,Rij,Self,I,Npop,FIT,ALLELES);
            vector<double> RChrom1_Fix = PopFix.getVec1();
            vector<double> RChrom2_Fix = PopFix.getVec4();

            int indPop(0);
            for(indPop = 0 ; indPop < Npop ; ++indPop)
            {
                if(RChrom1_Fix[indPop] != Rjk)
                {
                    NbAllMut += 1;
                }

                if(RChrom2_Fix[indPop] != Rjk)
                {
                    NbAllMut += 1;
                }
            }

            //cout << NbAllMut << endl;

        } while ((NbAllMut != 0) && (NbAllMut != 2*Npop));

        //cout << "Fixation done" << endl;

        if(NbAllMut == 2*Npop)
        {
            ComptFix += 1;
        }
    }

    if(Results)
    {
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "RESULTATS" << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Valeurs des paramètres : " << endl;
        Results << endl;
        Results << "Force initiale des Cis-acting factors : " << Ei << endl;
        Results << "Coefficient de dominance : " << h << endl;
        Results << "Intensité de la sélection au Locus A : " << S << endl;
        Results << "Taux de mutation au Locus A : " << mutA << endl;
        Results << "Taux de mutation au Locus E : " << mutE << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus E : " << SigE << endl;
        Results << "Ecart-type de la loi normale régissant les mutations multiplicatives au Locus R : " << SigR << endl;
        Results << "Taux de recombinaison entre les locus T et E : " << Rij << endl;
        Results << "Taux de recombinaison initial entre les locus E et A : " << Rjk << endl;
        Results << "Taux d'autofécondation : " << Self << endl;
        Results << "Intensité de la contrainte évolutive : " << I << endl;
        Results << "Taille de la population d'individus diploïdes : " << Npop << endl;
        Results << "Nombre de fixations simulées : " << Nit << endl;
        Results << "Nombre de générations simulées pour désinitialiser la population : " << TpsDesinit << endl;
        if(FIT == "Selection")
        {
            Results << "Sélection active" << endl;
        }
        else if(FIT == "Derive")
        {
            Results << "Dérive" << endl;
        }
        if(ALLELES == "2")
        {
            Results << "2 allèles au Locus A" << endl;
        }
        else if(ALLELES == "Infinite")
        {
            Results << "Nombre infini d'allèles au locus A" << endl;
        }
        Results << endl;
        Results << "_____________________________________________________________________________________________________________" << endl;
        Results << endl;
        Results << "Probabilité de fixation = " << ComptFix << endl;
        Results << endl;
    }
}
