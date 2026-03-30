//
// Created by admin on 2026/2/5.
//

#ifndef INC_1_NPA_H
#define INC_1_NPA_H
#include <lambda_lanczos/lambda_lanczos.hpp>
using lambda_lanczos::LambdaLanczos;

using namespace std;

template <typename T>
int kroneckerDelta(const T& a, const T& b);

void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove);

void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove);

void permuteWithSTL(vector<int>& nums, vector<vector<int>>& permutations);

int generateBasis(const vector<Pair*>& pairs, const int& pairNumber, const vector<int>& basis,
    vector<vector<int>>& bases, const vector<PairLimit*>& limits);

bool judgeBasis(const vector<PairLimit*>& limits, const vector<int>& basis);

int generateJi(const vector<Pair*>& pairs, const vector<vector<int>>& bases,
    vector<vector<vector<int>>>& JisList);

int generateBasesJ(const vector<vector<vector<int>>>& JisList, const vector<vector<int>>& bases,
    vector<vector<pair<int, int>>>& basesJ);

int generateJiForOne(const vector<Pair*>& pairs, const vector<int>& basis, const vector<int>& jis, const int& pairNumber,
    vector<vector<int>>& jisList);

int buildBlockJ(const vector<vector<pair<int, int>>>& basesJ, vector<pair<int, int>>& blockJ);

int initializeOrbitM(const vector<Orbit*>& orbits, vector<OrbitM*>& orbitMs);

int initializePairM(const vector<Pair*>& pairs, const vector<Orbit*>& orbits,
    const vector<OrbitM*>& orbitMs, vector<PairM*>& pairMs, vector<vector<int>>& pairJMMap);

int generateMBases(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<vector<int>>& bases, const vector<vector<int>>& pairJMMap, vector<vector<int>>& basesM0,
    vector<vector<int>>& basesM1);

int generateMBasis(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<int>& basis, const vector<int>& basisM, const vector<vector<int>>& pairJMMap,
    vector<vector<int>>& basesM0, vector<vector<int>>& basesM1);

int judgeMBasis(const vector<PairM*>& pairMs, const vector<int>& basisM);

int generalMultiCommutator(const vector<PairNew*>& bra, const vector<PairNew*>& ket,
    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>>& resultQ);

double overlapMSchemeOne(const vector<PairM*>& pairMs, const vector<int>& basisMBra,
    const vector<int>& basisMKet, const int orbitNumber);

int N2Qaby6(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar);

int N1Bab(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar, Eigen::SparseMatrix<double, Eigen::RowMajor>& Bab);

int overlapMScheme(const vector<PairM*>& pairMs, const vector<vector<int>>& basesM, const int orbitNumber,
    Eigen::MatrixXd& overlapMap);

int gramSchmidt(Eigen::MatrixXd& overlapMapJ,
    const vector<pair<int, int>>& blocks,
    vector<vector<pair<int, int>>>& basesJ,
    Eigen::MatrixXd& orthogonalBasis,
    Eigen::MatrixXd& transformMatrix0,
    Eigen::MatrixXd& transformMatrix1);

int matrixElementMtoJ(const Eigen::MatrixXd& overlapMapM,
    const Eigen::MatrixXd& transformMatrixBra,
    const Eigen::MatrixXd& transformMatrixKet,
    Eigen::MatrixXd& overlapMapJ);

int removeUselessBasesJ(vector<vector<pair<int, int>>>& basesJ,
    Eigen::MatrixXd& transformMatrix0,
    Eigen::MatrixXd& transformMatrix1);

int transferMatrix(const vector<Pair*>& pairs, const vector<PairM*>& pairMs,
    const vector<vector<pair<int, int>>>& basesJ, const vector<vector<int>>& basisM,
    Eigen::MatrixXd& transformMatrix);

double transferMatrixJMOne(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const vector<int>& basis,
    const vector<int>& jis, const vector<int>& basisM);

double CgJ0JnrM0m1mn(const int J0, const int M0, const int r, const int start,
    const vector<int>& jis, vector<int>& mis);

int gaby6(const vector<PairM*>& pairMs, const vector<int>& bra, const vector<int>& ket, const int orbitNumber,
    Eigen::MatrixXd& gaby6Matrix);

int fab(const vector<PairM*>& pairMs, const vector<int>& bra, const vector<int>& ket, const int orbitNumber,
    Eigen::MatrixXd& fabMatrix);

double reducedMatrixElement(const int& J, const int& JPrime, const int& M, const int& MPrime, const int& t,
    const int& u, const double& matrixElement);

int qab(const int orbitNumber, const vector<OrbitM*>& orbitMs, const vector<Orbit*>& orbits,
    Eigen::MatrixXd& qabMatrix);

int oaby6(const int orbitNumber, const vector<Orbit*>& orbitAC, const vector<OrbitM*>& orbitMAC,
    const vector<Orbit*>& orbitBD, const vector<OrbitM*>& orbitMBD,
    const map<TBMEJ, map<pair<int, int>, double>>& tbmeJMap, Eigen::SparseMatrix<double>& oaby6Matrix);

int indexOrbitJM(const vector<Orbit*>& orbits, const OrbitM* om);

int separateTBMEMap(map<TBMEJ, map<pair<int, int>, double>>& total, map<TBMEJ, map<pair<int, int>, double>>& ppMap,
    map<TBMEJ, map<pair<int, int>, double>>& pnMap, map<TBMEJ, map<pair<int, int>, double>>& nnMap,
    const int orbitNumberP);

int lanczos(Eigen::MatrixXd& matrix, vector<double>& eigenValues, vector<vector<double>>& eigenVectors);

int interactionAntiSymmetric(map<TBMEJ, map<pair<int, int>, double>>& imp, const vector<Orbit*>& orbits);

int calHamiltonianMatrix(const vector<PairM*>& pairMs, const vector<vector<int>>& basesM, const int orbitNumber,
    const vector<vector<Eigen::MatrixXd>>& fabMap,
    Eigen::SparseMatrix<double>& oaby6Matrix,
    Eigen::MatrixXd& qabMatrix,
    Eigen::MatrixXd& hamMatrix);

int separateReducedJ(const Eigen::MatrixXd& reducedJ, const vector<vector<pair<int, int>>>& basesJ,
    vector<pair<int, int>>& blocks, vector<pair<int, Eigen::MatrixXd>>& SRJBlocks);

int qabtJtoM(const vector<OrbitM*>& orbitMs,
    const vector<Orbit*>& orbits,
    const vector<int>& ts,
    const int& rank,
    const vector<vector<Eigen::MatrixXd>>& qPN,  // after svd, q^i(j_alpha j_beta t)
    vector<vector<Eigen::MatrixXd>>& qabtMMatrixMap);

int generateQMMatrix(const vector<vector<int>>& braBasesM,
    const vector<vector<int>>& ketBasesM,
    const vector<vector<Eigen::MatrixXd>>& qabtMMatrixMap,
    const vector<vector<Eigen::MatrixXd>>& fabMap,
    vector<vector<Eigen::MatrixXd>>& qMMatrix);

int generateFabMap(const vector<PairM*>& pairMs,
    const vector<vector<int>>& basesMBra,
    const vector<vector<int>>& basesMKet,
    const int orbitNumber,
    vector<vector<Eigen::MatrixXd>>& fabMap);

double QPiDotQNuJM(const vector<int>& bra, const vector<int>& ket, // 0 Pi, 1 Nu, 2 J, 3 M
    const vector<int>& ts, const int k, const int r, // r is i
    const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<Eigen::MatrixXd>>& QNuReducedElement,
    const vector<vector<Eigen::MatrixXd>>& QPiReducedElement);

void decomposeRightTensor(const Eigen::MatrixXd& mat,
                          const vector<int>& dims,
                          int rank,
                          std::vector<Eigen::MatrixXd>& q_pi,
                          Eigen::MatrixXd& V_it,
                          std::vector<Eigen::MatrixXd>& q_nu);

void getSVDResult(const map<TBMEJ, map<pair<int, int>, double>>& tbmePN,
                          const vector<Orbit*> orbitsP,
                          const vector<Orbit*> orbitsN,
                          const vector<int>& dims,
                          const vector<int>& ts,
                          vector<vector<Eigen::MatrixXd>>& q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
                          vector<Eigen::MatrixXd>& V_it,               // 输出: V_it(i, t)
                          vector<vector<Eigen::MatrixXd>>& q_nu,
                          const vector<double>& strengthVec);

int generateBasesJPNOne(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const int J,
    vector<vector<int>>& basesPN);

int generateTwoBodyPN(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& bra,
    const vector<vector<int>>& ket,
    const Eigen::MatrixXd& hamJP,
    const Eigen::MatrixXd& hamJN,
    Eigen::MatrixXd& HPiNuMatrixElement);

int generateTwoBodyPN(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& bra,
    const vector<vector<int>>& ket,
    const vector<int>& ts,
    const vector<Eigen::MatrixXd>& V_it,
    const vector<vector<Eigen::MatrixXd>>& qNuReducedElement,
    const vector<vector<Eigen::MatrixXd>>& qPiReducedElement,
    Eigen::MatrixXd& HPiNuMatrixElement);

bool isTriangle(int a, int b, int c);

int interactionAntiSymmetricPN(map<TBMEJ, map<pair<int, int>, double>>& imp);

template <typename T>
int initializeTwoDimVectors(int rows, int cols, vector<vector<T>>& twoDVector) {
    twoDVector.resize(rows);
    for (int i = 0; i < rows; i++) {
        twoDVector[i].resize(cols);
    }
    return 0;
}

int generateReducedOrthogonalQJ(const vector<int>& ts,
    const int rank,
    const vector<vector<pair<int, int>>>& basesJ,
    const vector<vector<Eigen::MatrixXd>>& qJMatrix00,
    const vector<vector<Eigen::MatrixXd>>& qJMatrix01,
    const Eigen::MatrixXd& orthogonalBasis,
    vector<vector<Eigen::MatrixXd>>& reducedJMatrix);


#endif //INC_1_NPA_H