#ifndef MODELRFIX_H_INCLUDED
#define MODELRFIX_H_INCLUDED
#include <string>

using namespace std;

void ModelRFix(double Ei, double mutE, double mutA, double SigR, double SigE, double h, double S, double Rij, double Rjk, double Self, double I, int Npop, int Nit, int TpsDesinit, string FIT, string ALLELES, string SUIVIFIT, string OutputPathway);

#endif // MODELRFIX_H_INCLUDED
