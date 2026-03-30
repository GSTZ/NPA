#include "CG.h"
#include "input.h"
#include "npa.h"

#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <ranges>
#include <set>
#include <chrono>

//using namespace std;

template <typename T>
int kroneckerDelta(const T& a, const T& b) {
    return a == b ? 1 : 0;
}

template <typename T, int Options>
T sparseTrace(const Eigen::SparseMatrix<T, Options>& mat) {
    T trace = T(0);
    for (Eigen::Index i = 0; i < min(mat.rows(), mat.cols()); ++i) {
        trace += mat.coeff(i, i);
    }
    return trace;
}


void removeRow(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    // 如果删除的不是最后一行，则将下面的行整体上移覆盖
    if (rowToRemove < numRows) {
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols)
            = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
    }

    // 缩小矩阵尺寸，但保留原有数据
    matrix.conservativeResize(numRows, numCols);
}


void removeColumn(Eigen::MatrixXd& matrix, unsigned int colToRemove) {
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

    // 如果删除的不是最后一列，则将右侧的列整体左移覆盖
    if (colToRemove < numCols) {
        matrix.block(0, colToRemove, numRows, numCols - colToRemove)
            = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
    }

    // 缩小矩阵尺寸，但保留原有数据
    matrix.conservativeResize(numRows, numCols);
}


void permuteWithSTL(vector<int>& nums, vector<vector<int>>& permutations) {
    sort(nums.begin(), nums.end()); // 必须先排序
    int count = 0;
    do {
        vector<int> tem;
        for (int x : nums) {
            tem.push_back(x);
        }
        permutations.push_back(tem);
        count++;
    } while (next_permutation(nums.begin(), nums.end()));
}



int generateBasis(const vector<Pair*>& pairs, const int& pairNumber, const vector<int>& basis,
    vector<vector<int>>& bases, const vector<PairLimit*>& limits) {
    if (basis.size() == pairNumber) {
        if (judgeBasis(limits, basis)) {
            bases.push_back(basis);
        }
    } else {
        int last, pairTypes, i;
        pairTypes = pairs.size();
        if (basis.empty()) last = 0;
        else last = basis.back();
        for (i = last; i < pairTypes; i++) {
            vector<int> basisNext = basis;
            basisNext.push_back(i);
            generateBasis(pairs, pairNumber, basisNext, bases, limits);
        }
    }
    return 0;
}


bool judgeBasis(const vector<PairLimit*>& limits, const vector<int>& basis) {
    bool flag = true;
    for (const auto& limit : limits) {
        int count = 0;
        for (const auto& index : limit->indexes) {
            count += std::count(basis.begin(), basis.end(), index);
        }
        if (count < limit->limitMin || count > limit->limitMax) {
            flag = false;
        }
    }
    return flag;
}


int generateJi(const vector<Pair*>& pairs, const vector<vector<int>>& bases, vector<vector<vector<int>>>& JisList) {
    for (int i = 0; i < bases.size(); i++) {
        auto basis = bases[i];
        vector<vector<int>> jisList;
        vector<int> jis;
        jis.push_back(pairs[basis[0]]->j);
        generateJiForOne(pairs, basis, jis, static_cast<int>(basis.size()), jisList);
        JisList.push_back(jisList);
    }
    return 0;
}


int generateBasesJ(const vector<vector<vector<int>>>& JisList, const vector<vector<int>>& bases,
    vector<vector<pair<int, int>>>& basesJ) {
    set<int> Js;
    for (const auto& jis : JisList) {
        for (const auto& basis : jis) {
            Js.insert(basis.back());
        }
    }

    for (const auto J : Js) {
        for (int i = 0; i < bases.size(); i++) {
            for (int j = 0; j < JisList[i].size(); j++) {
                if (JisList[i][j].back() == J) {
                    vector<pair<int, int>> temVector;
                    for (int k = 0; k < bases[i].size(); k++) {
                        temVector.emplace_back(bases[i][k], JisList[i][j][k]);
                    }
                    basesJ.push_back(temVector);
                }
            }
        }
    }
    return 0;
}


int generateJiForOne(const vector<Pair*>& pairs, const vector<int>& basis, const vector<int>& jis, const int& pairNumber,
    vector<vector<int>>& jisList) {
    if (jis.size() == pairNumber) {
        jisList.push_back(jis);
    } else {
        auto pairNext = pairs[basis[jis.size()]];
        auto ji = jis.back();
        vector<int> jiTem;
        if (ji == pairNext->j) {
            for (int i = abs(ji - pairNext->j); i <= abs(ji + pairNext->j); i += 2) {
                jiTem.push_back(i);
            }
        } else {
            for (int i = abs(ji - pairNext->j); i <= abs(ji + pairNext->j); i += 2) {
                jiTem.push_back(i);
            }
        }

        for (const auto& jt : jiTem) {
            auto jisNext = jis;
            jisNext.push_back(jt);
            generateJiForOne(pairs, basis, jisNext, pairNumber, jisList);
        }

    }
    return 0;
}


int buildBlockJ(const vector<vector<pair<int, int>>>& basesJ, vector<pair<int, int>>& blockJ) {
    blockJ.clear();
    int J0 = basesJ[0].back().second;
    int start = 0;
    for (int i = 0; i < basesJ.size(); i++) {
        if (basesJ[i].back().second != J0) {
            pair temPair = {start, i};
            blockJ.push_back(temPair);
            start = i;
            J0 = basesJ[i].back().second;
        }
    }
    pair temPair = {start, basesJ.size()};
    blockJ.emplace_back(temPair);
    return 0;
}


int initializeOrbitM(const vector<Orbit*>& orbits, vector<OrbitM*>& orbitMs) {
    for (const auto& orbit : orbits) {
        for (int i = -orbit->j; i <= orbit->j; i += 2) {
            auto om = new OrbitM(orbit->n, orbit->l, orbit->j, i, orbit->spe);
            orbitMs.push_back(om);
        }
    }
    return 0;
}


int initializePairM(const vector<Pair*>& pairs, const vector<Orbit*>& orbits,
    const vector<OrbitM*>& orbitMs, vector<PairM*>& pairMs, vector<vector<int>>& pairJMMap) {
    int index = 0;
    for (const auto& pair : pairs) {
        vector<int> ms;
        for (int i = -pair->j; i <= pair->j; i += 2) {
            Eigen::SparseMatrix<double, Eigen::RowMajor> pab(orbitMs.size(), orbitMs.size());
            for (int a = 0; a < orbitMs.size(); a++) {
                int ja;
                auto oma = orbitMs[a];
                for (ja = 0; ja < orbits.size(); ja++) {
                    if (orbits[ja]->n == oma->n && orbits[ja]->l == oma->l && orbits[ja]->j == oma->j) break;
                }
                for (int b = 0; b < orbitMs.size(); b++) {
                    auto omb = orbitMs[b];
                    int jb;
                    for (jb = 0; jb < orbits.size(); jb++) {
                        if (orbits[jb]->n == omb->n && orbits[jb]->l == omb->l && orbits[jb]->j == omb->j) break;
                    }
                    double v = pair->yabr.coeff(ja, jb) * CgInt(oma->j, oma->m, omb->j,
                        omb->m, pair->j, i);
                    if (abs(v) > 1e-15) {
                        pab.coeffRef(a, b) = v;
                    }
                }
            }
            auto pm = new PairM(index, pair->index, pair->j, i, pair->parity, pab);
            ms.push_back(index);
            pairMs.push_back(pm);
            ++index;
        }
        pairJMMap.push_back(ms);
    }
    return 0;
}

int generateMBases(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<vector<int>>& bases, const vector<vector<int>>& pairJMMap, vector<vector<int>>& basesM0,
    vector<vector<int>>& basesM1) {
    for (const auto& basis : bases) {
        vector<int> basisM;
        generateMBasis(pairs, pairMs, pairNumber, basis, basisM, pairJMMap, basesM0, basesM1);
    }
    return 0;
}

int generateMBasis(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<int>& basis, const vector<int>& basisM, const vector<vector<int>>& pairJMMap,
    vector<vector<int>>& basesM0, vector<vector<int>>& basesM1) {
    int index = basisM.size();
    if (index == pairNumber) {
        auto M = judgeMBasis(pairMs, basisM);
        if (M == 0) basesM0.push_back(basisM);
        else if (M == 2) basesM1.push_back(basisM);
    } else {
        auto pair = pairs[basis[index]];
        for (const auto pmIndex : pairJMMap[pair->index]) {
            auto pairM = pairMs[pmIndex];
            vector<int> basisMNext = basisM;
            basisMNext.push_back(pairM->index);
            generateMBasis(pairs, pairMs, pairNumber, basis, basisMNext, pairJMMap, basesM0, basesM1);
        }
        /*for (int i = 0; i < pairMs.size(); i++) {
            PairM pairM = pairMs[i];
            if (ranges::contains(pairJMMap[pair.index], pairM.index)) {
                vector<int> basisMNext = basisM;
                basisMNext.push_back(pairM.index);
                generateMBasis(pairs, pairMs, pairNumber, basis, basisMNext, pairJMMap, basesM);
            }
        }*/
    }
    return 0;
}


int judgeMBasis(const vector<PairM*>& pairMs, const vector<int>& basisM) {
    int M = 0;
    vector<int> order;
    for (const auto& bmi : basisM) {
        M += pairMs[bmi]->m;
        order.push_back(pairMs[bmi]->index);
    }

    if (!ranges::is_sorted(order)) {
        M = -1;
    }

    return M;
}


int overlapMScheme(const vector<PairM*>& pairMs, const vector<vector<int>>& basesM, const int orbitNumber,
    Eigen::MatrixXd& overlapMap) {
    overlapMap.resize(basesM.size(), basesM.size());
    overlapMap.setZero();
    for (int i = 0; i < basesM.size(); i++) {
        for (int j = 0; j < basesM.size(); j++) {
            //cout << i << endl;
            const auto& basisMBra = basesM[i];
            const auto& basisMKet = basesM[j];
            double overlap = overlapMSchemeOne(pairMs, basisMBra, basisMKet, orbitNumber);
            overlapMap(i, j) = overlap;
        }
    }


    /*auto basisMBra = basesM[0];
    auto basisMKet = basesM[2];
    double overlap = overlapMSchemeOne(pairMs, basisMBra, basisMKet, orbitNumber);
    cout << overlap << endl;

    basisMBra = basesM[2];
    basisMKet = basesM[0];
    overlap = overlapMSchemeOne(pairMs, basisMBra, basisMKet, orbitNumber);
    cout << overlap << endl;*/

    return 0;
}

double overlapMSchemeOne(const vector<PairM*>& pairMs, const vector<int>& basisMBra,
    const vector<int>& basisMKet, const int orbitNumber) {
    double overlap = 0.0;
    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>> resultQ;
    vector<PairNew*> bra, ket;
    // bra and ket have the same length.
    for (int i = 0; i < basisMBra.size() - 1; i++) {
        auto pNewBra = new PairNew();
        pNewBra->pab = pairMs[basisMBra[i]]->pab;
        bra.push_back(pNewBra);
    }

    for (int i = 0; i < basisMKet.size(); i++) {
        auto pNewKet = new PairNew();
        pNewKet->pab = pairMs[basisMKet[i]]->pab;
        ket.push_back(pNewKet);
    }
    //generalMultiCommutator(bra, ket, resultQ);

    //vector<vector<Eigen::SparseMatrix<double, Eigen::RowMajor>>> qbar;
    Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber * orbitNumber, orbitNumber * orbitNumber);

    //N2Qaby6(bra, ket, orbitNumber, qbar);


    /*for (int i = 0; i < orbitNumber; i++) {
        for (int j = 0; j < orbitNumber; j++) {
            for (int k = 0; k < orbitNumber; k++) {
                for (int l = 0; l < orbitNumber; l++) {
                    for (int m = 0; m < resultQ.size(); m++) {
                        qbar.coeffRef(i * orbitNumber + k, j * orbitNumber + l) +=
                            (resultQ[m].lhs().coeff(i, j) * resultQ[m].rhs().coeff(k, l)
                            - resultQ[m].lhs().coeff(i, k) * resultQ[m].rhs().coeff(j, l)
                            + resultQ[m].lhs().coeff(i, l) * resultQ[m].rhs().coeff(j, k)
                            + resultQ[m].lhs().coeff(j, k) * resultQ[m].rhs().coeff(i, l)
                            - resultQ[m].lhs().coeff(j, l) * resultQ[m].rhs().coeff(i, k)
                            + resultQ[m].lhs().coeff(j, l) * resultQ[m].rhs().coeff(i, k)) / 6;
                    }
                }
            }
        }
    }*/

    Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);
    //auto pN1 = pairMs[basisMBra[basisMBra.size() - 2]]->pab;

    N1Bab(bra, ket, orbitNumber, qbar, Bab);

    /*for (int i = 0; i < orbitNumber; i++) {
        for (int j = 0; j < orbitNumber; j++) {
            for (int k = 0; k < orbitNumber; k++) {
                for (int l = 0; l < orbitNumber; l++) {
                    Bab.coeffRef(i * orbitNumber + k, j * orbitNumber + l) += 12 * pN1.coeff(i, j)
                    * qbar.coeff(i * orbitNumber + k, j * orbitNumber + l);
                }
            }
        }
    }*/

    Eigen::SparseMatrix<double, Eigen::RowMajor> pnBn = pairMs[basisMBra.back()]->pab * Bab;
    overlap = -2 * sparseTrace(pnBn);

    return overlap;
}


int N2Qaby6(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar) {

    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>> resultQ;

    generalMultiCommutator(bra, ket, resultQ);


    for (int i = 0; i < orbitNumber; i++) {
        for (int j = 0; j < orbitNumber; j++) {
            for (int k = 0; k < orbitNumber; k++) {
                for (int l = 0; l < orbitNumber; l++) {
                    for (int m = 0; m < resultQ.size(); m++) {
                        qbar.coeffRef(i * orbitNumber + k, j * orbitNumber + l) +=
                            (resultQ[m].lhs().coeff(i, j) * resultQ[m].rhs().coeff(k, l)
                            - resultQ[m].lhs().coeff(i, k) * resultQ[m].rhs().coeff(j, l)
                            + resultQ[m].lhs().coeff(i, l) * resultQ[m].rhs().coeff(j, k)
                            + resultQ[m].lhs().coeff(j, k) * resultQ[m].rhs().coeff(i, l)
                            - resultQ[m].lhs().coeff(j, l) * resultQ[m].rhs().coeff(i, k)
                            + resultQ[m].lhs().coeff(k, l) * resultQ[m].rhs().coeff(i, j)) / 6;
                    }
                }
            }
        }
    }
    return 0;
}


int N1Bab(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar, Eigen::SparseMatrix<double, Eigen::RowMajor>& Bab) {

    auto pN1 = bra.back()->pab;
    auto braNext = bra;
    braNext.pop_back();

    N2Qaby6(braNext, ket, orbitNumber, qbar);

    for (int i = 0; i < orbitNumber; i++) {
        for (int j = 0; j < orbitNumber; j++) {
            for (int k = 0; k < orbitNumber; k++) {
                for (int l = 0; l < orbitNumber; l++) {
                    Bab.coeffRef(i, j) += 12 * pN1.coeff(k, l)
                    * qbar.coeff(i * orbitNumber + k, j * orbitNumber + l);
                }
            }
        }
    }

    return 0;
}


int generalMultiCommutator(const vector<PairNew*>& bra, const vector<PairNew*>& ket,
    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>>& resultQ) {
    auto N = ket.size();
    // recursion end
    if (bra.empty() && N == 2) {
        auto Q = ket[0]->pab * ket[1]->pab;
        resultQ.emplace_back(Q);
    } else {
        for (int k = 0; k < N; k++) {
            if (k != N - 1) { // contract the trace to P_{N}
                auto braNext = bra;
                auto pnPrime = braNext.back();
                braNext.pop_back();
                auto ketNext = ket;
                auto pn = ketNext.back();
                auto pk = ketNext[k];
                ketNext.pop_back();
                ketNext.erase(ketNext.begin() + k);
                auto Pk = new PairNew();
                Eigen::SparseMatrix<double> pTem = pk->pab * pnPrime->pab;
                Pk->pab = -2 * sparseTrace(pTem) * pn->pab;
                ketNext.push_back(Pk);
                generalMultiCommutator(braNext, ketNext, resultQ);
            } else { // contract the trace to P_{N-1} when k = N
                auto braNext = bra;
                auto pnPrime = braNext.back();
                braNext.pop_back();
                auto ketNext = ket;
                auto pnk = ketNext.back();
                ketNext.pop_back();
                auto pnMinus1 = ketNext.back();
                ketNext.pop_back();
                auto Pk = new PairNew();
                Eigen::SparseMatrix<double> pTem = pnk->pab * pnPrime->pab;
                Pk->pab = -2 * sparseTrace(pTem) * pnMinus1->pab;
                ketNext.push_back(Pk);
                generalMultiCommutator(braNext, ketNext, resultQ);
            }
        }
        for (int k = 0; k < N; k++) {
            for (int i = 0; i < k; i++) {
                auto braNext = bra;
                auto pnPrime = braNext.back();
                braNext.pop_back();
                auto ketNext = ket;
                auto pk = ketNext[k];
                auto pi = ketNext[i];
                auto Pik = new PairNew();
                // i < k, so erase order can not be changed.
                ketNext.erase(ketNext.begin() + k);
                ketNext.erase(ketNext.begin() + i);
                Pik->pab = 4 * (pk->pab * pnPrime->pab * pi->pab + pi->pab * pnPrime->pab * pk->pab);
                ketNext.push_back(Pik);
                generalMultiCommutator(braNext, ketNext, resultQ);
            }
        }
    }
    return 0;
}


int gramSchmidt(Eigen::MatrixXd& overlapMapJ,
    const vector<pair<int, int>>& blocks,
    vector<vector<pair<int, int>>>& basesJ,
    Eigen::MatrixXd& orthogonalBasis,
    Eigen::MatrixXd& transformMatrix0,
    Eigen::MatrixXd& transformMatrix1) {
    int dim = static_cast<int>(overlapMapJ.cols());

    Eigen::MatrixXd omNew;
    omNew.resize(dim, dim);
    omNew.setZero();
    for (int i = 0; i < dim; i++) {
        omNew(i, i) = 1.0;
    }

    vector<int> validIndices;
    vector<int> invalidIndices;

    // 预计算每个向量属于哪个块
    vector block_index(dim, -1);
    for (int b = 0; b < blocks.size(); b++) {
        int start = blocks[b].first;
        int end = blocks[b].second;
        for (int i = start; i < end; i++) {  // 注意：i < end
            block_index[i] = b;
        }
    }

    // 对每个块单独处理
    for (const auto& block : blocks) {
        int start = block.first;
        int end = block.second;

        for (int k = start; k < end; k++) {  // 注意：k < end
            // 投影减法
            for (int j : validIndices) {
                // 只处理同一块内的向量
                if (block_index[j] != block_index[k]) continue;

                double overlap_kj = 0.0;
                // 只计算块内的重叠积分
                for (int alpha = start; alpha < end; alpha++) {  // alpha < end
                    for (int beta = start; beta < end; beta++) {  // beta < end
                        overlap_kj += omNew(k, alpha) * overlapMapJ(alpha, beta) * omNew(j, beta);
                    }
                }

                // 向量减法
                for (int alpha = 0; alpha < dim; alpha++) {
                    omNew(k, alpha) -= overlap_kj * omNew(j, alpha);
                }
            }

            // 计算范数（只在块内计算）
            double norm_k_sq = 0.0;
            for (int alpha = start; alpha < end; alpha++) {  // alpha < end
                for (int beta = start; beta < end; beta++) {  // beta < end
                    norm_k_sq += omNew(k, alpha) * overlapMapJ(alpha, beta) * omNew(k, beta);
                }
            }

            if (norm_k_sq < 1e-10) {
                //cerr << "Warning: Found near-zero vector at k=" << k << ", possibly due to linear dependence.\n";
                invalidIndices.push_back(k);
                continue;
            }

            double norm_k = sqrt(norm_k_sq);
            for (int alpha = 0; alpha < dim; alpha++) {
                omNew(k, alpha) /= norm_k;
            }
            validIndices.push_back(k);
        }
    }

    // 构建正交基
    //vector<vector<double>> orthogonalBasis;
    //for (int idx : validIndices) {
        //orthogonalBasis.push_back(omNew[idx]);
    //}
    orthogonalBasis.resize(validIndices.size(), omNew.cols());
    for (Eigen::Index i = 0; i < validIndices.size(); ++i) {
        orthogonalBasis.row(i) = omNew.row(validIndices[i]);
    }

    //orthogonalBasis = omNew(validIndices, Eigen::all);

    // 删除无效基矢（从后往前删除以避免索引问题）
    ranges::sort(invalidIndices, greater<int>());
    for (int k : invalidIndices) {
        basesJ.erase(basesJ.begin() + k);
        //ystrallp.erase(ystrallp.begin() + k);
        removeRow(overlapMapJ, k);
        removeColumn(overlapMapJ, k);

        removeColumn(orthogonalBasis, k);

        removeRow(transformMatrix0, k);
        removeRow(transformMatrix1, k);
    }

    return 0;
}


int matrixElementMtoJ(const Eigen::MatrixXd& overlapMapM,
    const Eigen::MatrixXd& transformMatrixBra,
    const Eigen::MatrixXd& transformMatrixKet,
    Eigen::MatrixXd& overlapMapJ) {
    overlapMapJ = transformMatrixBra * overlapMapM * transformMatrixKet.transpose();
    return 0;
}


int removeUselessBasesJ(vector<vector<pair<int, int>>>& basesJ,
    Eigen::MatrixXd& transformMatrix0,
    Eigen::MatrixXd& transformMatrix1) {
    for (int i = transformMatrix0.rows() - 1; i > -1; i--) {
        if (transformMatrix0.row(i).isZero(1e-13)) {
            removeRow(transformMatrix0, i);
            removeRow(transformMatrix1, i);
            basesJ.erase(basesJ.begin() + i);
        }
    }
    return 0;
};


int transferMatrix(const vector<Pair*>& pairs, const vector<PairM*>& pairMs,
    const vector<vector<pair<int, int>>>& basesJ, const vector<vector<int>>& basisM,
    Eigen::MatrixXd& transformMatrix) {
    int rows = basesJ.size();
    const int cols = static_cast<int>(basisM.size());

    transformMatrix.resize(rows, cols);
    transformMatrix.setZero();

    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            vector<int> basisPair;
            vector<int> jis;
            const auto& basis = basesJ[j];
            for (int k = 0; k < basis.size(); k++) {
                basisPair.push_back(basis[k].first);
                jis.push_back(basis[k].second);
            }
            transformMatrix(j, i) = transferMatrixJMOne(pairs, pairMs, basisPair,
                jis, basisM[i]);
        }
    }

    return 0;
}


double transferMatrixJMOne(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const vector<int>& basis,
    const vector<int>& jis, const vector<int>& basisM) {

    for (int i = 0; i < basis.size(); i++) {
        if (pairs[basis[i]]->index != pairMs[basisM[i]]->indexJ) {
            return 0.0;
        }
    }

    double Tcg = 1.0;
    vector<pair<int, int>> classify;
    int type = pairs[basis[0]]->index;
    pair se = {0, 0};
    for (const int bs : basis) {
        if (pairs[bs]->index == type) {
            se.second += 1;
        } else {
            classify.push_back(se);
            se = {se.second, se.second + 1};
            type = pairs[bs]->index;
        }
    }
    classify.push_back(se);

    int M0 = 0;
    for (int i = 0; i < classify.size(); i++) {
        auto cf = classify[i];
        vector jisSliced(jis.begin() + cf.first, jis.begin() + cf.second);
        vector<int> mis;
        for (int j = cf.first; j < cf.second; j++) {
            mis.push_back(pairMs[basisM[j]]->m);
        }

        int J0;
        if (i == 0) {
            J0 = 0;
        } else {
            J0 = jis[cf.first - 1];
        }

        Tcg *= CgJ0JnrM0m1mn(J0, M0, pairs[basis[cf.first]]->j, cf.first, jisSliced, mis);

        M0 += accumulate(mis.begin(), mis.end(), 0);
    }

    return Tcg;
}


double CgJ0JnrM0m1mn(const int J0, const int M0, const int r, const int start,
    const vector<int>& jis, vector<int>& mis) {
    double cg = 0.0;
    vector<vector<int>> permutation;
    permuteWithSTL(mis, permutation);
    for (const auto& element : permutation) {
        double cgTem = 1.0;
        auto J = J0;
        auto M = M0;
        for (int i = 0; i < element.size(); i++) {
            cgTem *= CgInt(J, M, r, element[i], jis[i], M + element[i]);
            M += element[i];
            J = jis[i];
        }
        cg += cgTem;
    }
    return cg;
}


int gaby6(const vector<PairM*>& pairMs, const vector<int>& bra, const vector<int>& ket, const int orbitNumber,
    Eigen::MatrixXd& gaby6Matrix) {
    gaby6Matrix.resize(orbitNumber * orbitNumber, orbitNumber * orbitNumber);
    gaby6Matrix.setZero();
    vector<PairNew*> ketNext;
    for (int i = 0; i < ket.size(); i++) {
        auto pw = new PairNew();
        pw->pab = pairMs[ket[i]]->pab;
        ketNext.push_back(pw);
    }

    for (int i = 0; i < bra.size(); i++) {
        auto braTem = bra;
        braTem.erase(braTem.begin() + i);
        auto Pk = pairMs[bra[i]];
        vector<PairNew*> braNext;
        for (int j = 0; j < braTem.size(); j++) {
            auto pw = new PairNew();
            pw->pab = pairMs[braTem[j]]->pab;
            braNext.push_back(pw);
        }
        Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber * orbitNumber, orbitNumber * orbitNumber);
        Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

        N1Bab(braNext, ketNext, orbitNumber, qbar, Bab);
        for (int alpha = 0; alpha < orbitNumber; alpha++) {
            for (int beta = 0; beta < orbitNumber; beta++) {
                auto pkab = Pk->pab.coeff(alpha, beta);
                for (int gamma = 0; gamma < orbitNumber; gamma++) {
                    for (int delta = 0; delta < orbitNumber; delta++) {
                        gaby6Matrix(alpha * orbitNumber + gamma, beta * orbitNumber + delta)
                        += 4 * pkab * Bab.coeff(gamma, delta);
                    }
                }
            }
        }
    }

    for (int k = 0; k < bra.size(); k++) {
        for (int i = 0; i < k; i++) {
            auto Pi = pairMs[bra[i]];
            auto Pk = pairMs[bra[k]];

            auto braTem = bra;
            braTem.erase(braTem.begin() + k);
            braTem.erase(braTem.begin() + i);

            vector<PairNew*> braNext;
            for (int j = 0; j < braTem.size(); j++) {
                auto pw = new PairNew();
                pw->pab = pairMs[braTem[j]]->pab;
                braNext.push_back(pw);
            }

            Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber * orbitNumber,
                                                              orbitNumber * orbitNumber);
            N2Qaby6(braNext, ketNext, orbitNumber, qbar);

            for (int gamma = 0; gamma < orbitNumber; gamma++) {
                for (int delta = 0; delta < orbitNumber; delta++) {
                    Eigen::SparseMatrix<double, Eigen::RowMajor> qbarAB(orbitNumber, orbitNumber);

                    for (int alpha = 0; alpha < orbitNumber; alpha++) {
                        for (int beta = 0; beta < orbitNumber; beta++) {
                            double val = qbar.coeff(alpha * orbitNumber + gamma,
                                                    beta * orbitNumber + delta);
                            if (val != 0.0) {
                                qbarAB.insert(alpha, beta) = val;
                            }
                        }
                    }

                    Eigen::SparseMatrix<double, Eigen::RowMajor> Sab = Pk->pab * qbarAB * Pi->pab
                    + Pi->pab * qbarAB * Pk->pab;
                    Sab *= 96.0;

                    for (int alpha = 0; alpha < orbitNumber; alpha++) {
                        for (int beta = 0; beta < orbitNumber; beta++) {
                            gaby6Matrix(alpha * orbitNumber + gamma, beta * orbitNumber + delta)
                                += Sab.coeff(alpha, beta);
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int fab(const vector<PairM*>& pairMs, const vector<int>& bra, const vector<int>& ket, const int orbitNumber,
    Eigen::MatrixXd& fabMatrix) {
    fabMatrix.resize(orbitNumber, orbitNumber);
    fabMatrix.setZero();
    vector<PairNew*> ketNext;
    for (int i = 0; i < ket.size(); i++) {
        auto pw = new PairNew();
        pw->pab = pairMs[ket[i]]->pab;
        ketNext.push_back(pw);
    }

    for (int i = 0; i < bra.size(); i++) {
        auto braTem = bra;
        braTem.erase(braTem.begin() + i);
        auto Pk = pairMs[bra[i]];
        vector<PairNew*> braNext;
        for (int j = 0; j < braTem.size(); j++) {
            auto pw = new PairNew();
            pw->pab = pairMs[braTem[j]]->pab;
            braNext.push_back(pw);
        }
        Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber * orbitNumber, orbitNumber * orbitNumber);
        Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

        N1Bab(braNext, ketNext, orbitNumber, qbar, Bab);
        fabMatrix += 4 * Pk->pab * Bab.transpose();
    }
    return 0;
}


double reducedMatrixElement(const int& J, const int& JPrime, const int& M, const int& MPrime, const int& t,
    const int& u, const double& matrixElement) {
    double result = matrixElement / CgInt(JPrime, MPrime, t, u, J, M);
    return result;
}


int qab(const int orbitNumber, const vector<OrbitM*>& orbitMs, const vector<Orbit*>& orbits,
    Eigen::MatrixXd& qabMatrix) {
    qabMatrix.resize(orbitNumber, orbitNumber);
    qabMatrix.setZero();
    for (int i = 0; i < orbitNumber; ++i) {
        for (int j = 0; j < orbitNumber; ++j) {
            if (indexOrbitJM(orbits, orbitMs[i]) == indexOrbitJM(orbits, orbitMs[j]))
            qabMatrix(i, j) = orbitMs[i]->spe;
        }
    }
    return 0;
}


int oaby6(const int orbitNumber, const vector<Orbit*>& orbitAC, const vector<OrbitM*>& orbitMAC,
    const vector<Orbit*>& orbitBD, const vector<OrbitM*>& orbitMBD,
    const map<TBMEJ, map<pair<int, int>, double>>& tbmeJMap, Eigen::SparseMatrix<double>& oaby6Matrix) {
    oaby6Matrix.resize(orbitNumber * orbitNumber, orbitNumber * orbitNumber);
    oaby6Matrix.setZero();
    for (int alpha = 0; alpha < orbitNumber; alpha++) {
        for (int beta = 0; beta < orbitNumber; beta++) {
            for (int gamma = 0; gamma < orbitNumber; gamma++) {
                for (int delta = 0; delta < orbitNumber; delta++) {
                    const auto ma = orbitMAC[alpha]->m;
                    const auto mb = orbitMBD[beta]->m;
                    const auto mc = orbitMAC[gamma]->m;
                    const auto md = orbitMBD[delta]->m;
                    if (ma + mb == mc + md) {
                        double value = 0.0;
                        auto a = indexOrbitJM(orbitAC, orbitMAC[alpha]);
                        auto b = indexOrbitJM(orbitBD, orbitMBD[beta]);
                        auto c = indexOrbitJM(orbitAC, orbitMAC[gamma]);
                        auto d = indexOrbitJM(orbitBD, orbitMBD[delta]);
                        auto ja = orbitMAC[alpha]->j;
                        auto jb = orbitMBD[beta]->j;
                        auto jc = orbitMAC[gamma]->j;
                        auto jd = orbitMBD[delta]->j;
                        auto tj = TBMEJ(a, b, c, d);
                        if (tbmeJMap.find(tj) != tbmeJMap.end()) {
                            auto JTList = tbmeJMap.find(tj)->second;
                            for (const auto& jt : JTList) {
                                value += sqrt((1 + kroneckerDelta(a, b)) * (1 + kroneckerDelta(c, d))) * jt.second
                                * CgInt(ja, ma, jb, mb, jt.first.first, ma + mb)
                                * CgInt(jc, mc, jd, md, jt.first.first, ma + mb);
                            }
                            oaby6Matrix.coeffRef(alpha * orbitNumber + gamma, beta * orbitNumber + delta) += value;
                        }
                    }
                }
            }
        }
    }
    return 0;
}


int indexOrbitJM(const vector<Orbit*>& orbits, const OrbitM* om) {
    int index = -1;
    for (int i = 0; i < orbits.size(); i++) {
        auto os = orbits[i];
        if (os->n == om->n && os->l == om->l && os->j == om->j) {
            index = i;
            break;
        }
    }
    return index;
}


int separateTBMEMap(map<TBMEJ, map<pair<int, int>, double>>& total, map<TBMEJ, map<pair<int, int>, double>>& ppMap,
    map<TBMEJ, map<pair<int, int>, double>>& pnMap, map<TBMEJ, map<pair<int, int>, double>>& nnMap,
    const int orbitNumberP) {
    for (const auto& tt : total) {
        if (tt.first.a > orbitNumberP) {
            nnMap.emplace(tt.first, tt.second);
        } else if (tt.first.a <= orbitNumberP && tt.first.b > orbitNumberP) {
            pnMap.emplace(tt.first, tt.second);
        } else {
            ppMap.emplace(tt.first, tt.second);
        }
    }
    return 0;
}


int lanczos(Eigen::MatrixXd& matrix, vector<double>& eigenValues, vector<vector<double>>& eigenVectors) {
    auto mv_mul = [&](const vector<double>& in, vector<double>& out) {
        auto eigen_in = Eigen::Map<const Eigen::VectorXd>(&in[0], in.size());
        auto eigen_out = Eigen::Map<Eigen::VectorXd>(&out[0], out.size());

        eigen_out = matrix * eigen_in; // Easy version
    };

    LambdaLanczos<double> engine(mv_mul, matrix.cols(), false, 2);
    engine.run(eigenValues, eigenVectors);
    return 0;
}


int interactionAntiSymmetric(map<TBMEJ, map<pair<int, int>, double>>& imp, const vector<Orbit*>& orbits) {
    map<TBMEJ, map<pair<int, int>, double>> result;
    for (auto& ip : imp) {
        auto old = ip.first;
        int j1 = orbits[old.a]->j;
        int j2 = orbits[old.b]->j;
        int j3 = orbits[old.c]->j;
        int j4 = orbits[old.d]->j;
        result[old] = ip.second;

        TBMEJ tj = TBMEJ{old.c, old.d, old.a, old.b};
        result[tj] = ip.second;

        tj = TBMEJ{old.b, old.a, old.c, old.d};
        map<pair<int, int>, double> newJTV;
        for (auto& jtv : ip.second) {
            int J = jtv.first.first;
            int T = jtv.first.second;
            double v = jtv.second;
            double phase = -1 * pow(-1,(j1 + j2 - J) / 2);
            newJTV[jtv.first] = phase * v;
        }
        result[tj] = newJTV;

        tj = TBMEJ{old.c, old.d, old.b, old.a};
        result[tj] = newJTV;

        tj = TBMEJ{old.a, old.b, old.d, old.c};
        newJTV.clear();
        for (auto& jtv : ip.second) {
            int J = jtv.first.first;
            int T = jtv.first.second;
            double v = jtv.second;
            double phase = -1 * pow(-1,(j3 + j4 - J) / 2);
            newJTV[jtv.first] = phase * v;
        }
        result[tj] = newJTV;

        tj = TBMEJ{old.d, old.c, old.a, old.b};
        result[tj] = newJTV;

        tj = TBMEJ{old.b, old.a, old.d, old.c};
        newJTV.clear();
        for (auto& jtv : ip.second) {
            int J = jtv.first.first;
            int T = jtv.first.second;
            double v = jtv.second;
            double phase = pow(-1,(j1 + j2 + j3 + j4) / 2);
            newJTV[jtv.first] = phase * v;
        }
        result[tj] = newJTV;

        tj = TBMEJ{old.d, old.c, old.b, old.a};
        result[tj] = newJTV;

    }

    imp = result;

    return 0;
}

int interactionAntiSymmetricPN(map<TBMEJ, map<pair<int, int>, double>>& imp) {
    map<TBMEJ, map<pair<int, int>, double>> result;
    for (auto& ip : imp) {
        auto old = ip.first;
        result[old] = ip.second;

        TBMEJ tj = TBMEJ{old.c, old.d, old.a, old.b};
        result[tj] = ip.second;
    }

    imp = result;

    return 0;
}


int calHamiltonianMatrix(const vector<PairM*>& pairMs, const vector<vector<int>>& basesM, const int orbitNumber,
    const vector<vector<Eigen::MatrixXd>>& fabMap,
    Eigen::SparseMatrix<double>& oaby6Matrix,
    Eigen::MatrixXd& qabMatrix,
    Eigen::MatrixXd& hamMatrix) {
    hamMatrix.resize(basesM.size(), basesM.size());
    hamMatrix.setZero();
    Eigen::MatrixXd gaby6Matrix;
    for (int i = 0; i < basesM.size(); i++) {
        for (int j = 0; j < basesM.size(); j++) {
            const auto& basisMBra = basesM[i];
            const auto& basisMKet = basesM[j];
            gaby6(pairMs, basisMBra, basisMKet, orbitNumber, gaby6Matrix);
            Eigen::MatrixXd fabMatrix = fabMap[i][j];
            //fab(pairMs, basisMBra, basisMKet, orbitNumber, fabMatrix);
            double element = 0.0;
            for (int k = 0; k < gaby6Matrix.rows(); k++) {
                for (int l = 0; l < gaby6Matrix.cols(); l++) {
                    element += oaby6Matrix.coeffRef(k, l) * gaby6Matrix(k, l) / 4;
                }
            }
            for (int k = 0; k < qabMatrix.rows(); k++) {
                for (int l = 0; l < qabMatrix.cols(); l++) {
                    element += fabMatrix.coeffRef(k, l) * qabMatrix(k, l);
                }
            }
            hamMatrix(i, j) = element;
        }
    }
    return 0;
}


int separateReducedJ(const Eigen::MatrixXd& reducedJ, const vector<vector<pair<int, int>>>& basesJ,
    vector<pair<int, int>>& blocks, vector<pair<int, Eigen::MatrixXd>>& SRJBlocks) {
    for (const auto& b : blocks) {
        int J = basesJ[b.first].back().second;
        Eigen::MatrixXd oneBlock;
        int length = b.second - b.first;
        oneBlock.resize(length, length);
        oneBlock.setZero();
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                oneBlock(i, j) = reducedJ(i + b.first, j + b.first);
            }
        }
        pair p = {J, oneBlock};
        SRJBlocks.push_back(p);
    }
    return 0;
}


int generateHPiMu(const vector<int>& ts,
                          const vector<int>& dims,
                          const vector<vector<Eigen::MatrixXd>>& q_pi,
                          const vector<Eigen::MatrixXd>& V_it,
                          const vector<vector<Eigen::MatrixXd>>& q_nu,
                          Eigen::MatrixXd& hPiMu) {
    const int d0 = dims[0]; // j_alpha
    const int d1 = dims[1]; // j_beta
    const int d2 = dims[2]; // j_gamma
    const int d3 = dims[3]; // j_delta
    hPiMu.resize(d0 * d1, d2 * d3);
    hPiMu.setZero();
    for (int i = 0; i < ts.size(); i++) {
        int t = ts[i];
        const vector<Eigen::MatrixXd>& q_pi_temp = q_pi[i];
        const Eigen::MatrixXd& V_it_temp = V_it[t];
        const vector<Eigen::MatrixXd>& q_nu_temp = q_nu[i];
        for (int j = 0; j < q_pi_temp.size(); j++) {
            for (int alpha = 0; alpha < d0; alpha++) {
                for (int beta = 0; beta < d1; beta++) {
                    for (int gamma = 0; gamma < d2; gamma++) {
                        for (int delta = 0; delta < d3; delta++) {
                            hPiMu(alpha * d1 + beta, gamma * d3 + delta) += std::pow(-1, t / 2.0)
                            * sqrt(t + 1) * q_pi_temp[j](alpha, beta) * q_nu_temp[j](gamma, delta)
                            * V_it_temp(j, j);
                        }
                    }
                }
            }
        }
    }
    return 0;
}



void getSVDResult(const map<TBMEJ, map<pair<int, int>, double>>& tbmePN,
                          const vector<Orbit*> orbitsP,
                          const vector<Orbit*> orbitsN,
                          const vector<int>& dims,
                          const vector<int>& ts,
                          vector<vector<Eigen::MatrixXd>>& q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
                          vector<Eigen::MatrixXd>& V_it,               // 输出: V_it(i, t)
                          vector<vector<Eigen::MatrixXd>>& q_nu,
                          const vector<double>& strengthVec) {

    const int d0 = dims[0]; // j_alpha
    const int d1 = dims[1]; // j_beta
    const int d2 = dims[2]; // j_gamma
    const int d3 = dims[3]; // j_delta

    for (auto t : ts) {
        Eigen::MatrixXd mat(d0 * d1, d2 * d3);
        mat.setZero();
        for (const auto& tpn : tbmePN) {
            auto tj = tpn.first;
            auto alpha = tj.a;
            auto gamma = tj.b;
            auto beta = tj.c;
            auto delta = tj.d;
            auto ja = orbitsP[alpha]->j;
            auto jb = orbitsP[beta]->j;
            auto jc = orbitsN[gamma]->j;
            auto jd = orbitsN[delta]->j;
            for (const auto& jtv : tpn.second) {
                auto J = jtv.first.first;
                mat(alpha * d1 + beta, gamma * d3 + delta) += std::pow(-1, (J + jc + jb) / 2.0)
                * (J + 1) * sixjSymbols(ja / 2.0, jc / 2.0, J / 2.0, jd / 2.0, jb / 2.0, t / 2.0)
                * jtv.second;
            }
        }

        // 2. 使用 SVD 分解
        vector<Eigen::MatrixXd> q_pi_temp;
        Eigen::MatrixXd V_it_temp;
        vector<Eigen::MatrixXd> q_nu_temp;

        // 假设秩为10
        decomposeRightTensor(mat, dims, 100, q_pi_temp, V_it_temp, q_nu_temp);

        // 3. 存储结果
        q_pi.push_back(q_pi_temp);
        V_it.push_back(V_it_temp * strengthVec[0]);
        q_nu.push_back(q_nu_temp);
    }

}


// 假设四维张量的维度为 (n, n, n, n)
void decomposeRightTensor(const Eigen::MatrixXd& mat,
                          const vector<int>& dims,
                          int rank,
                          std::vector<Eigen::MatrixXd>& q_pi,
                          Eigen::MatrixXd& V_it,
                          std::vector<Eigen::MatrixXd>& q_nu) {

    const int d0 = dims[0]; // j_alpha
    const int d1 = dims[1]; // j_beta
    const int d2 = dims[2]; // j_gamma
    const int d3 = dims[3]; // j_delta

    // 2. SVD分解
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // 3. 确定有效秩
    const int max_rank = std::min(mat.rows(), mat.cols());
    if (rank > max_rank) {
        std::cerr << "Warning: Reducing rank from " << rank
                  << " to max possible " << max_rank << std::endl;
        rank = max_rank;
    }

    // 4. 提取成分
    Eigen::MatrixXd U = svd.matrixU().leftCols(rank);
    Eigen::MatrixXd V = svd.matrixV().leftCols(rank);
    Eigen::VectorXd S = svd.singularValues().head(rank);


    // 5. 重组因子矩阵
    q_pi.resize(rank);
    q_nu.resize(rank);
    for (int i = 0; i < rank; ++i) {
        // 左因子重组为 d0 x d1
        q_pi[i] = Eigen::Map<Eigen::MatrixXd>(U.col(i).data(), d0, d1);
        // q_pi[i]= q_pi[i].transpose();
        q_pi[i].transposeInPlace();

        // 右因子重组为 d2 x d3
        q_nu[i] = Eigen::Map<Eigen::MatrixXd>(V.col(i).data(), d2, d3);
        q_nu[i].transposeInPlace();
    }

    // 6. 奇异值矩阵
    V_it = S.asDiagonal();
}


int generateBasesJPNOne(const vector<vector<pair<int, int> > > &basesJP,
                        const vector<vector<pair<int, int> > > &basesJN,
                        const int J,
                        vector<vector<int> > &basesPN) {
    for (int i = 0; i < basesJN.size(); ++i) {
        for (int j = 0; j < basesJP.size(); ++j) {
            auto JPi = basesJP[j].back().second;
            auto JNu = basesJN[i].back().second;
            if (isTriangle(JPi, JNu, J)) {
                vector<int> one;
                one.push_back(j);
                one.push_back(i);
                one.push_back(J);
                basesPN.push_back(one);
            }
        }
    }
    return 0;
}


/*int generateTwoBodyPN(const vector<vector<pair<int, int>>>& BasesJP,
    const vector<vector<pair<int, int>>>& BasesJN,
    const vector<pair<int, int>>& bra,
    const vector<pair<int, int>>& ket,
    const vector<int>& ts,
    const vector<int>& dims,
    vector<vector<Eigen::MatrixXd>>& q_pi,  // 输出: q_pi[i](j_alpha, j_beta)
    vector<Eigen::MatrixXd>& V_it,               // 输出: V_it(i, t)
    vector<vector<Eigen::MatrixXd>>& q_nu) {

    const int d0 = dims[0]; // j_alpha
    const int d1 = dims[1]; // j_beta
    const int d2 = dims[2]; // j_gamma
    const int d3 = dims[3]; // j_delta

    for (int i = 0; i < bra.size(); ++i) {
        for (int j = 0; j < ket.size(); ++j) {
            auto braMuJ = BasesJN[bra[i].first];
            auto braPiJ = BasesJP[bra[i].second];
            auto ketMuJ = BasesJN[ket[j].first];
            auto ketPiJ = BasesJP[ket[j].second];
            Eigen::MatrixXd mat(d0 * d1, d2 * d3);
            mat.setZero();
            for (int k = 0; k < ts.size(); ++k) {
                auto t = ts[k];
                auto qPiTem = q_pi[k];
                auto qNuTem = q_nu[k];
                auto rank = V_it[k].size();
                for (int l = 0; l < rank; ++l) {
                    Eigen::MatrixXd mat1(d0 * d1, d2 * d3);
                    mat.setZero();
                    auto qabPi = qPiTem[l];
                    auto qabNu = qNuTem[l];
                }
            }
        }
    }


    auto rank = q_pi[0].size();
    for (int i = 0; i < rank; ++i) {
        for (int k = 0; k < ts.size(); ++k) {
            auto t = ts[k];
            for (int i1 = 0; i1 < bra.size(); ++i1) {
                for (int i2 = 0; i2 < ket.size(); ++i2) {

                }
            }
        }
    }

    return 0;
}*/

int generateTwoBodyPN(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& bra,
    const vector<vector<int>>& ket,
    const Eigen::MatrixXd& hamJP,
    const Eigen::MatrixXd& hamJN,
    Eigen::MatrixXd& HPiNuMatrixElement) {
    for (int i = 0; i < bra.size(); ++i) {
        for (int j = 0; j < ket.size(); ++j) {
            const auto& braNext = bra[i];
            const auto& ketNext = ket[j];
            int braJPIndex = braNext[0];
            int ketJPIndex = ketNext[0];
            int braJNIndex = braNext[1];
            int ketJNIndex = ketNext[1];
            if (basesJN[braJNIndex].back().second == basesJN[ketJNIndex].back().second && basesJP[braJPIndex].back().second == basesJP[ketJPIndex].back().second) {
                if (braJPIndex == ketJPIndex) HPiNuMatrixElement(i, j) += hamJN(braJNIndex, ketJNIndex);
                if (braJNIndex == ketJNIndex) HPiNuMatrixElement(i, j) += hamJP(braJPIndex, ketJPIndex);
            }
        }
    }
    return 0;
}


int generateTwoBodyPN(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& bra,
    const vector<vector<int>>& ket,
    const vector<int>& ts,
    const vector<Eigen::MatrixXd>& V_it,
    const vector<vector<Eigen::MatrixXd>>& qNuReducedElement,
    const vector<vector<Eigen::MatrixXd>>& qPiReducedElement,
    Eigen::MatrixXd& HPiNuMatrixElement) {

    HPiNuMatrixElement.resize(bra.size(), ket.size());
    HPiNuMatrixElement.setZero();

    auto rank = V_it[0].rows();
    for (int i = 0; i < rank; ++i) {
        for (int k = 0; k < ts.size(); ++k) {
            auto t = ts[k];
            for (int i1 = 0; i1 < bra.size(); ++i1) {
                for (int i2 = 0; i2 < ket.size(); ++i2) {
                    const auto& braNext = bra[i1];
                    const auto& ketNext = ket[i2];
                    HPiNuMatrixElement(i1, i2) += QPiDotQNuJM(braNext, ketNext, ts, k, i, basesJP,
                        basesJN, qNuReducedElement, qPiReducedElement)
                        * std::pow(-1, t / 2.0) * V_it[k](i, i);
                }
            }
        }
    }

    return 0;
}


int qabtJtoM(const vector<OrbitM*>& orbitMs,
    const vector<Orbit*>& orbits,
    const vector<int>& ts,
    const int& rank,
    const vector<vector<Eigen::MatrixXd>>& qPN,  // after svd, q^i(j_alpha j_beta t)
    vector<vector<Eigen::MatrixXd>>& qabtMMatrixMap) {
    auto orbitNumber = orbitMs.size();
    for (int r = 0; r < rank; ++r) {
        for (int k = 0; k < ts.size(); ++k) {
            auto t = ts[k];
            const auto& qPnTem = qPN[k][r];
            Eigen::MatrixXd qmMatrix;
            qmMatrix.resize(orbitNumber, orbitNumber);
            qmMatrix.setZero();
            for (int i = 0; i < orbitNumber; ++i) {
                auto bra = orbitMs[i];
                auto JAlpha = indexOrbitJM(orbits, bra);
                for (int j = 0; j < orbitNumber; ++j) {
                    auto ket = orbitMs[j];
                    auto JBeta = indexOrbitJM(orbits, ket);
                    qmMatrix(i, j) = std::pow(-1, (ket->j + ket->m) / 2.0)
                    * CgInt(bra->j, bra->m, ket->j, -ket->m, t, bra->m - ket->m)
                    * qPnTem(JAlpha, JBeta);
                }
            }
            qabtMMatrixMap[k][r] = qmMatrix;
        }
    }
    return 0;
}


int generateQMMatrix(const vector<vector<int>>& braBasesM,
    const vector<vector<int>>& ketBasesM,
    const vector<vector<Eigen::MatrixXd>>& qabtMMatrixMap,
    const vector<vector<Eigen::MatrixXd>>& fabMap,
    vector<vector<Eigen::MatrixXd>>& qMMatrix) {

    for (int m = 0; m < qabtMMatrixMap.size(); ++m) {
        for (int n = 0; n < qabtMMatrixMap[m].size(); ++n) {
            auto& oneBlock = qMMatrix[m][n];
            oneBlock.resize(braBasesM.size(), ketBasesM.size());
            oneBlock.setZero();
            auto qabMatrix = qabtMMatrixMap[m][n];
            for (int i = 0; i < braBasesM.size(); i++) {
                for (int j = 0; j < ketBasesM.size(); j++) {
                    const Eigen::MatrixXd& fabMatrix = fabMap[i][j];
                    double element = 0.0;
                    for (int k = 0; k < qabMatrix.rows(); k++) {
                        for (int l = 0; l < qabMatrix.cols(); l++) {
                            element += fabMatrix(k, l) * qabMatrix(k, l);
                        }
                    }
                    oneBlock(i, j) = element;
                }
            }
        }
    }
    return 0;
}



int generateFabMap(const vector<PairM*>& pairMs,
    const vector<vector<int>>& basesMBra,
    const vector<vector<int>>& basesMKet,
    const int orbitNumber,
    vector<vector<Eigen::MatrixXd>>& fabMap) {
    fabMap.resize(basesMBra.size());
    for (int i = 0; i < basesMBra.size(); ++i) {
        fabMap[i].resize(basesMKet.size());
    }
    Eigen::MatrixXd fabMatrix;
    for (int i = 0; i < basesMBra.size(); i++) {
        for (int j = 0; j < basesMKet.size(); j++) {
            const auto& basisMBra = basesMBra[i];
            const auto& basisMKet = basesMKet[j];
            fab(pairMs, basisMBra, basisMKet, orbitNumber, fabMatrix);
            fabMap[i][j] = fabMatrix;
        }
    }
    return 0;
}


double QPiDotQNuJM(const vector<int>& bra, const vector<int>& ket, // 0 Pi, 1 Nu, 2 J, 3 M
    const vector<int>& ts, const int k, const int r, // r is i
    const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<Eigen::MatrixXd>>& QNuReducedElement,
    const vector<vector<Eigen::MatrixXd>>& QPiReducedElement) {
    double result = 0.0;
    int JPiPrime = basesJP[bra[0]].back().second;
    int JNuPrime = basesJN[bra[1]].back().second;
    int JPi = basesJP[ket[0]].back().second;
    int JNu = basesJN[ket[1]].back().second;
    int J = bra[2];
    int t = ts[k];
    result = std::pow(-1, (JPiPrime + JNu + J + t) / 2.0)
           * sixjSymbols(JNuPrime / 2.0, JPiPrime / 2.0, J / 2.0, JPi / 2.0, JNu / 2.0, t / 2.0)
           * sqrt(JPiPrime + 1) * sqrt(JNuPrime + 1)
           * QPiReducedElement[k][r](bra[0], ket[0])
           * QNuReducedElement[k][r](bra[1], ket[1]);
    return result;
}


bool isTriangle(int a, int b, int c) {
    return a >= 0 && b >= 0 && c >= 0 &&
           a + b >= c &&
           a + c >= b &&
           b + c >= a;
}


int generateReducedOrthogonalQJ(const vector<int>& ts,
    const int rank,
    const vector<vector<pair<int, int>>>& basesJ,
    const vector<vector<Eigen::MatrixXd>>& qJMatrix00,
    const vector<vector<Eigen::MatrixXd>>& qJMatrix01,
    const Eigen::MatrixXd& orthogonalBasis,
    vector<vector<Eigen::MatrixXd>>& reducedJMatrix) {
    for (int i = 0; i < ts.size(); i++) {
        for (int j = 0; j < rank; j++) {
            auto t = ts[i];
            Eigen::MatrixXd oneBlock;
            oneBlock.resize(basesJ.size(), basesJ.size());
            oneBlock.setZero();
            for (int k = 0; k < basesJ.size(); k++) {
                for (int l = 0; l < basesJ.size(); l++) {
                    int J = basesJ[k].back().second;
                    int JPrime = basesJ[l].back().second;
                    if (isTriangle(JPrime, J, t)) {
                        if ((J + JPrime + t) % 4 == 0) {
                            oneBlock(k, l) = qJMatrix00[i][j](k, l) / CgInt(JPrime, 0, t, 0, J, 0);
                        } else {
                            oneBlock(k, l) = qJMatrix01[i][j](k, l) / CgInt(JPrime, 2, t, -2, J, 0);
                        }
                    }
                }
            }
            oneBlock = orthogonalBasis * oneBlock * orthogonalBasis.transpose();
            reducedJMatrix[i][j] = oneBlock;
        }
    }
    return 0;
}


