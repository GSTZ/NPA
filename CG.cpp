//
// Created by admin on 2026/2/3.
//
#include <cmath>

using namespace std;

int factorial(int n) {
    if (n == 0) return 1;
    int result = 1;
    for (int i=1; i<=n; i++) {
        result *= i;
    }
    return result;
}

double factorial_double(int n) {
    return factorial(n);
}

double Cg(const double j1, const double m1, const double j2, const double m2, const double J, const double M) {
    double CGs;
    if (j1 < 0 || j2 < 0 || J < 0) {
        CGs=0;
    }
    else if (j1 < m1 || m1 < -j1) {
        CGs=0;
    }
    else if (j2 < m2 || m2 < -j2) {
        CGs=0;
    }
    else if (J < M || M < -J) {
        CGs=0;
    }
    else if (j1 - m1 - trunc(j1 - m1) != 0) {
        CGs=0;
    }
    else if (j2 - m2 - trunc(j2-m2) != 0) {
        CGs=0;
    }
    else if (J - M - trunc(J-M) != 0) {
        CGs=0;
    }
    else if (m1 + m2 != M) {
        CGs=0;
    }
    else if ((j1 + j2 < J) || (J < abs(j1 - j2))) {
        CGs=0;
    }
    else {
        int f, fj1, fm1, fj2, fm2, fJ, fM;
        double b1, b2;
        fj1=static_cast<int>(j1 * 2);
        fm1=static_cast<int>(m1 * 2);
        fj2=static_cast<int>(j2 * 2);
        fm2=static_cast<int>(m2 * 2);
        fJ=static_cast<int>(J * 2);
        fM=static_cast<int>(M * 2);


        b1 = (fJ + 1) * factorial_double((fj1 + fj2 - fJ) / 2) / factorial_double((fj1 + fj2 + fJ + 2) / 2)
        * factorial_double((fj1 - fm1) / 2)/factorial_double((fj1 - fj2 + fJ) / 2)
        * factorial_double((fj2 - fm2) / 2)/factorial_double((fj2 - fj1 + fJ) / 2)
        * factorial_double((fJ + fM) / 2)/factorial_double((fj1 + fm1) / 2)
        * factorial_double((fJ - fM) / 2)/factorial_double((fj2 + fm2) / 2);
        b2=0.0;
        for (f = max(0,(fJ-fj2-fm1)/2); f <= (fJ+fj1+fj2)/2; ++f) {
            if (j2+J-m1-f>=0 && j1-m1-f>=0 && J-M-f>=0 && j2-J+m1+f>=0) {
                if ((f+(fj1-fm1)/2) % 2 == 0) {
                    b2 = b2 + factorial_double((fj1+fm1+2*f)/2)
                    /factorial_double((fj1-fm1-2*f)/2)*factorial_double((fj2+fJ-fm1-2*f)/2)
                    /factorial_double((fJ-fM-2*f)/2)/factorial_double((fj2-fJ+fm1+2*f)/2)
                    /factorial_double(f);
                }
                else {
                    b2 = b2-factorial_double((fj1+fm1+2*f)/2)
                    /factorial_double((fj1-fm1-2*f)/2)*factorial_double((fj2+fJ-fm1-2*f)/2)
                    /factorial_double((fJ-fM-2*f)/2)/factorial_double((fj2-fJ+fm1+2*f)/2)
                    /factorial_double(f);
                }
            }
        }
        CGs=sqrt(b1)*b2;
    }
    return CGs;
}


double threejSymbols(double j1, double j2, double j3, double m1, double m2, double m3) {
    double threeJ, m;
    m = -m3;
    if (static_cast<int>(j1-j2-m3) % 2 == 0) {
        threeJ = Cg(j1,m1,j2,m2,j3,m)/sqrt(2*j3+1);
    }
    else {
        threeJ = -1 * Cg(j1,m1,j2,m2,j3,m)/sqrt(2*j3+1);
    }
    return threeJ;
}



double sixjSymbols(double j1, double j2, double j3, double m1, double m2, double m3) {
    double sixJ, a1, a2, a3, a4, b1, b2;
    int fj1, fj2, fj3, fm1, fm2, fm3, fc1, fc2, fc3, fc4, fc5, fc6, fc7, f, fcmax, fcmin, ff;

    if (j1 < 0 || j2 < 0 || j3 < 0 || m1 < 0 || m2 < 0 || m3 < 0) sixJ = 0;
    else if (j1 + j2 - j3 - trunc(j1 + j2 - j3) != 0) sixJ = 0;
    else if (j1 + j2 < j3 || j3 < abs(j1 - j2)) sixJ = 0;
    else if (j2 + m1 - m3 - trunc(j2 + m1 - m3) != 0) sixJ = 0;
    else if (j2 + m1 < m3 || m3 < abs(m1-j2)) sixJ = 0;
    else if (j3 + m1 - m2 - trunc(j3 + m1 - m2) != 0) sixJ = 0;
    else if (j3 + m1 < m2 || m2 < abs(m1 - j3)) sixJ = 0;
    else if (j1 + m3 - m2 - trunc(j1 + m3 - m2) != 0) sixJ = 0;
    else if (j1 + m3 < m2 || m2 < abs(j1 - m3)) sixJ = 0;
    else {
        fj1 = static_cast<int>(j1*2);
        fm1 = static_cast<int>(m1*2);
        fj2 = static_cast<int>(j2*2);
        fm2 = static_cast<int>(m2*2);
        fj3 = static_cast<int>(j3*2);
        fm3 = static_cast<int>(m3*2);

        a1 = sqrt(factorial_double((fj1+fj2-fj3)/2)*factorial_double((fj2+fj3-fj1)/2)
        *factorial_double((fj3+fj1-fj2)/2)/factorial_double((fj1+fj2+fj3+2)/2));
        a2 = sqrt(factorial_double((fj1+fm2-fm3)/2)*factorial_double((fm2+fm3-fj1)/2)
        *factorial_double((fm3+fj1-fm2)/2)/factorial_double((fj1+fm2+fm3+2)/2));
        a3 = sqrt(factorial_double((fm1+fj2-fm3)/2)*factorial_double((fj2+fm3-fm1)/2)
        *factorial_double((fm3+fm1-fj2)/2)/factorial_double((fm1+fj2+fm3+2)/2));
        a4 = sqrt(factorial_double((fm1+fm2-fj3)/2)*factorial_double((fm2+fj3-fm1)/2)
        *factorial_double((fj3+fm1-fm2)/2)/factorial_double((fm1+fm2+fj3+2)/2));
        b1 = a1 * a2 * a3 * a4;

        fc1=fj1+fm2+fm3;
        fc2=fm1+fj2+fm3;
        fc3=fm1+fm2+fj3;
        fc4=fj1+fj2+fj3;
        fc5=fj1+fj2+fm1+fm2;
        fc6=fj3+fj2+fm3+fm2;
        fc7=fj1+fj3+fm1+fm3;
        fcmax = fc2 > fc1 ? fc2 : fc1;
        fcmax = fc3 > fcmax ? fc3 : fcmax;
        fcmax = fc4 > fcmax ? fc4 : fcmax;
        fcmin = fc5 < fc6 ? fc5 : fc6;
        fcmin = fc7 < fcmin ? fc7 : fcmin;
        f = fcmax % 2 == 0 ? fcmax / 2 : fcmax / 2 + 1;
        b2 = 0.0;
        while (f * 2 <= fcmin) {
            ff = 2 * f;
            if (f % 2 == 0) {
                b2 = b2 + factorial_double(f+1)/factorial_double((ff-fc1)/2)
                / factorial_double((ff-fc2)/2)/factorial_double((ff-fc3)/2)/factorial_double((ff-fc4)/2)
                / factorial_double((fc5-ff)/2)/factorial_double((fc6-ff)/2)/factorial_double((fc7-ff)/2);
            } else {
                b2 = b2 - factorial_double(f+1) / factorial_double((ff-fc1)/2)
                / factorial_double((ff-fc2)/2) / factorial_double((ff-fc3)/2) / factorial_double((ff-fc4)/2)
                / factorial_double((fc5-ff)/2) / factorial_double((fc6-ff)/2) / factorial_double((fc7-ff)/2);
            }
            f = ++f;
        }
        sixJ = b1 * b2;
    }
    return sixJ;
}

double CgInt(const int j1, const int m1, const int j2, const int m2, const int J, const int M) {
    return Cg(static_cast<double>(j1) / 2.0, static_cast<double>(m1) / 2.0,
        static_cast<double>(j2) / 2.0, static_cast<double>(m2) / 2.0,
        static_cast<double>(J) / 2.0, static_cast<double>(M) / 2.0);
}


