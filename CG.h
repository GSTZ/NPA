//
// Created by admin on 2026/2/3.
//

#ifndef INC_1_CG_H
#define INC_1_CG_H

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int factorial(double n);
double factorial_double(int n);
double double_factorial(double x);
double Cg(const double j1, const double m1, const double j2, const double m2, const double J, const double M);
double threejSymbols(double j1, double j2, double j3, double m1, double m2, double m3);
double sixjSymbols(double j1, double j2, double j3, double m1, double m2, double m3);
double CgInt(const int j1, const int m1, const int j2, const int m2, const int J, const int M);

double sixjSymbolsInt(int j1, int j2, int j3, int m1, int m2, int m3);

#endif //INC_1_CG_H