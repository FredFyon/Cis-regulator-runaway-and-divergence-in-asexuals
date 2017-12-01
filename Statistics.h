#ifndef STATISTICS_H
#define STATISTICS_H
#include <vector>

using namespace std;

class Statistic
{
    public:
    Statistic(vector<double> x1, vector<double> x2);
    double getMean() const;
    double getVariance() const;
    double getStandardDeviation() const;

    private:
    vector<double> m1;
    vector<double> m2;
};

#endif // STATISTICS_H
