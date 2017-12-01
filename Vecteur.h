#ifndef VECTEUR_H
#define VECTEUR_H

using namespace std;

class Vecteur
{
    public:
    Vecteur(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10);
    double getx1() const;
    double getx2() const;
    double getx3() const;
    double getx4() const;
    double getx5() const;
    double getx6() const;
    double getx7() const;
    double getx8() const;
    double getx9() const;
    double getx10() const;

    private:
    double m1;
    double m2;
    double m3;
    double m4;
    double m5;
    double m6;
    double m7;
    double m8;
    double m9;
    double m10;

};

#endif // VECTEUR_H
