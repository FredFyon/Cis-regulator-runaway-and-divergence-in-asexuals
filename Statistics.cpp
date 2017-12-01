#include <string>
#include <vector>
#include <math.h>
#include "Statistics.h"

using namespace std;

Statistic::Statistic(vector<double> x1, vector<double> x2):m1(x1),m2(x2)
{

}

double Statistic::getMean() const
{
    double size = m1.size();
    double moy(0);
    int i(0);
    for(i = 0 ; i < size ; i++)
    {
        moy += m1[i]+m2[i];
    }
    moy = moy/(2.*size);

    return moy;
}

double Statistic::getVariance() const
{
    double size = m1.size();
    Statistic V(m1,m2);
    double mean = V.getMean();
    double variance(0);
    int i(0);
    for(i = 0 ; i < size ; i++)
    {
        variance += ((m1[i] + m2[i])/2. - mean)*((m1[i] + m2[i])/2. - mean);
    }
    variance /= size;

    return variance;
}

double Statistic::getStandardDeviation() const
{
    Statistic V(m1,m2);
    double variance = V.getVariance();
    double standardDeviation = sqrt(variance);

    return standardDeviation;
}
