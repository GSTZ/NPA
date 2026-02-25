//
// Created by admin on 2026/2/5.
//

#ifndef INC_1_NPA_H
#define INC_1_NPA_H
#include <vector>

using namespace std;

void permuteWithSTL(vector<int>& nums, vector<vector<int>>& permutations);

int generateBasis(const vector<Pair*>& pairs, const int& pairNumber, const vector<int>& basis,
    vector<vector<int>>& bases, const vector<PairLimit*>& limits);

bool judgeBasis(const vector<PairLimit*>& limits, const vector<int>& basis);

int generateJi(const vector<Pair*>& pairs, const vector<vector<int>>& bases,
    vector<vector<vector<int>>>& JisList);

int generateJiForOne(const vector<Pair*>& pairs, const vector<int>& basis, const vector<int>& jis, const int& pairNumber,
    vector<vector<int>>& jisList);

int initializeOrbitM(const vector<Orbit*>& orbits, vector<OrbitM*>& orbitMs);

int initializePairM(const vector<Pair*>& pairs, const vector<Orbit*>& orbits,
    const vector<OrbitM*>& orbitMs, vector<PairM*>& pairMs, vector<vector<int>>& pairJMMap);

int generateMBases(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<vector<int>>& bases, const vector<vector<int>>& pairJMMap, vector<vector<int>>& basesM);

int generateMBasis(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<int>& basis, const vector<int>& basisM, const vector<vector<int>>& pairJMMap,
    vector<vector<int>>& basesM);

bool judgeMBasis(const vector<PairM*>& pairMs, const vector<int>& basisM);

int generalMultiCommutator(const vector<PairNew*>& bra, const vector<PairNew*>& ket,
    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>>& resultQ);

double overlapMSchemeOne(const vector<PairM*>& pairMs, const vector<int>& basisMBra,
    const vector<int>& basisMKet, const int orbitNumber);

int N2Qaby6(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar);

int N1Bab(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar, Eigen::SparseMatrix<double, Eigen::RowMajor>& pN1,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& Bab);

int overlapMScheme(const vector<PairM*>& pairMs, const vector<vector<int>>& basesM, const int orbitNumber,
    vector<vector<double>>& overlapMap);

double transferMatrixJMOne(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const vector<int>& basis,
    const vector<int>& jis, const vector<int>& basisM);

double CgJ0JnrM0m1mn(const int J0, const int M0, const int r, const vector<int>& jis, vector<int> mis);

#endif //INC_1_NPA_H