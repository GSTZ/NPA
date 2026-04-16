//
// Created by admin on 2026/2/3.
//
#include <cmath>
#include <iostream>
#include <algorithm>
#include "CGC.h"

using namespace std;

// 双阶乘函数 (支持分数)
double double_factorial(double x) {
    if (x <= 0) return 1; // 双阶乘停止条件
    double result = 1.0;
    for (double i = x; i > 0; i -= 2) {
        result *= i;
    }
    return result;
}

/*double CgInt(const int j1, const int m1, const int j2, const int m2, const int J, const int M) {
    return Cg(static_cast<double>(j1) / 2.0, static_cast<double>(m1) / 2.0,
        static_cast<double>(j2) / 2.0, static_cast<double>(m2) / 2.0,
        static_cast<double>(J) / 2.0, static_cast<double>(M) / 2.0);
}

double sixjSymbolsInt(int j1, int j2, int j3, int m1, int m2, int m3) {
    return sixjSymbols(j1 / 2.0, j2 / 2.0, j3 / 2.0, m1 / 2.0, m2 / 2.0, m3 / 2.0);
}*/


double CgInt(int j1, int m1, int j2, int m2, int J, int M) {
    auto cg2 = util::CG(j1, j2, J, m1, m2, M);
    return cg2;
}

double sixjSymbolsInt(int j1, int j2, int j3, int j4, int j5, int j6) {
    return util::wigner_6j(j1, j2, j3, j4, j5, j6);
}


