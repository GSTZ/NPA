#include "CG.h"
#include "input.h"
#include "npa.h"

#include <string>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <ranges>
#include <set>

using namespace std;

template <typename T, int Options>
T sparseTrace(const Eigen::SparseMatrix<T, Options>& mat) {
    T trace = T(0);
    for (Eigen::Index i = 0; i < std::min(mat.rows(), mat.cols()); ++i) {
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
            for (int i = abs(ji - pairNext->j); i <= abs(ji + pairNext->j); i += 4) {
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


int initializeOrbitM(const vector<Orbit*>& orbits, vector<OrbitM*>& orbitMs) {
    for (const auto& orbit : orbits) {
        for (int i = -orbit->j; i <= orbit->j; i += 2) {
            auto om = new OrbitM(orbit->n, orbit->l, orbit->j, i);
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
                    pab.insert(a, b) = pair->yabr.coeff(ja, jb) * CgInt(oma->j, oma->m, omb->j,
                        omb->m, pair->j, i);
                }
            }
            auto pm = new PairM(index, pair->j, i, pair->parity, pab);
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
    for (const auto& bmi : basisM) {
        M += pairMs[bmi]->m;
    }
    return M;
}


int overlapMScheme(const vector<PairM*>& pairMs, const vector<vector<int>>& basesM, const int orbitNumber,
    Eigen::MatrixXd& overlapMap) {

    for (int i = 0; i < basesM.size(); i++) {
        for (int j = 0; j < basesM.size(); j++) {
            const auto& basisMBra = basesM[i];
            const auto& basisMKet = basesM[j];
            double overlap = overlapMSchemeOne(pairMs, basisMBra, basisMKet, orbitNumber);
            overlapMap(i, j) = overlap;
        }
    }

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
    Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber ^ 2, orbitNumber ^ 2);

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

    for (int i = 0; i < basisMBra.size(); i++) {
        //delete bra[i];
        delete ket[i];
    }
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
                    Bab.coeffRef(i * orbitNumber + k, j * orbitNumber + l) += 12 * pN1.coeff(i, j)
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


int gramSchmidt(Eigen::MatrixXd& overlapMapJ, const vector<pair<int, int>>& blocks,
    vector<vector<pair<int, int>>>& basesJ, Eigen::MatrixXd& orthogonalBasis) {
    int dim = static_cast<int>(overlapMapJ.size());

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
                cerr << "Warning: Found near-zero vector at k=" << k << ", possibly due to linear dependence.\n";
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
    }

    return 0;
}


int matrixElementMtoJ(const Eigen::MatrixXd& overlapMapM, const Eigen::MatrixXd& transformMatrix,
    Eigen::MatrixXd& overlapMapJ) {
    overlapMapJ = transformMatrix * overlapMapM * transformMatrix.transpose();
    return 0;
}


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
        int J0 = jis[cf.first - 1];

        Tcg *= CgJ0JnrM0m1mn(J0, M0, pairs[basis[cf.first]]->j, jisSliced, mis);

        M0 += accumulate(mis.begin(), mis.end(), 0);
    }

    return Tcg;
}


double CgJ0JnrM0m1mn(const int J0, const int M0, const int r, const vector<int>& jis, vector<int> mis) {
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
        }
        cg += cgTem;
    }
    return cg;
}


int gaby6(const vector<PairM*>& pairMs, const vector<int>& bra, const vector<int>& ket, const int orbitNumber,
    Eigen::MatrixXd& gaby6Matrix) {
    gaby6Matrix.resize(orbitNumber ^ 2, orbitNumber ^ 2);
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
        Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber ^ 2, orbitNumber ^ 2);
        Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

        N1Bab(braNext, ketNext, orbitNumber, qbar, Bab);
        for (int alpha = 0; alpha < orbitNumber; alpha++) {
            for (int beta = 0; beta < orbitNumber; beta++) {
                for (int gamma = 0; gamma < orbitNumber; gamma++) {
                    for (int delta = 0; delta < orbitNumber; delta++) {
                        auto pkab = Pk->pab.coeff(alpha, beta);
                        gaby6Matrix(alpha * orbitNumber + gamma, beta * orbitNumber + delta)
                        += 4 * pkab * Bab.coeff(gamma, delta);
                    }
                }
            }
        }
    }

    for (int gamma = 0; gamma < orbitNumber; gamma++) {
        for (int delta = 0; delta < orbitNumber; delta++) {
            Eigen::SparseMatrix<double, Eigen::RowMajor> Sab(orbitNumber, orbitNumber);
            for (int k = 1; k < bra.size(); k++) {
                for (int i = 0; i < k; i++) {
                    auto braTem = bra;
                    auto Pi = pairMs[bra[i]];
                    auto Pk = pairMs[bra[k]];
                    braTem.erase(braTem.begin() + k);
                    braTem.erase(braTem.begin() + i);
                    vector<PairNew*> braNext;
                    for (int j = 0; j < braTem.size(); j++) {
                        auto pw = new PairNew();
                        pw->pab = pairMs[braTem[j]]->pab;
                        braNext.push_back(pw);
                    }
                    Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber ^ 2, orbitNumber ^ 2);
                    N2Qaby6(braNext, ketNext, orbitNumber, qbar);
                    Eigen::SparseMatrix<double, Eigen::RowMajor> qbarAB(orbitNumber, orbitNumber);
                    for (int alpha = 0; alpha < orbitNumber; alpha++) {
                        for (int beta = 0; beta < orbitNumber; beta++) {
                            if (qbar.coeff(alpha * orbitNumber + gamma, beta * orbitNumber + delta) != 0.0) {
                                qbarAB.insert(alpha, beta) = qbar.coeff(alpha * orbitNumber + gamma,
                                    beta * orbitNumber + delta);
                            }
                        }
                    }
                    Sab += Pk->pab * qbarAB * Pi->pab + Pi->pab * qbarAB * Pk->pab;
                }
            }
            Sab *= 96.0;
            for (int alpha = 0; alpha < orbitNumber; alpha++) {
                for (int beta = 0; beta < orbitNumber; beta++) {
                    gaby6Matrix(alpha * orbitNumber + gamma, beta * orbitNumber + delta)
                    += Sab.coeff(alpha, beta);
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
        Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber ^ 2, orbitNumber ^ 2);
        Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

        N1Bab(braNext, ketNext, orbitNumber, qbar, Bab);
        fabMatrix += Pk->pab * Bab.transpose();
    }
    return 0;
}


double reducedMatrixElement(const int& J, const int& JPrime, const int& M, const int& MPrime, const int& t,
    const int& u, const double& matrixElement) {
    double result = matrixElement / CgInt(JPrime, MPrime, t, u, J, M);
    return result;
}


int qab(const int orbitNumber, const vector<OrbitM*>& orbitMs, Eigen::MatrixXd& qabMatrix) {
    qabMatrix.resize(orbitNumber, orbitNumber);
    qabMatrix.setZero();
    for (int i = 0; i < orbitNumber; ++i) {
        qabMatrix(i, i) = orbitMs[i]->spe;
    }
    return 0;
}


int oaby6(const int orbitNumber, const vector<Orbit*>& orbitAC, const vector<OrbitM*>& orbitMAC,
    const vector<Orbit*>& orbitBD, const vector<OrbitM*>& orbitMBD,
    const map<TBMEJ, map<pair<int, int>, double>>& tbmeJMap, Eigen::SparseMatrix<double>& oaby6Matrix) {
    oaby6Matrix.resize(orbitNumber ^ 2, orbitNumber ^ 2);
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
                        auto JTList = tbmeJMap.find(tj)->second;
                        for (const auto& jt : JTList) {
                            value += sqrt((1 + kroneckerDelta(a, b)) * (1 + kroneckerDelta(c, d))) * jt.second
                            * CgInt(ja, ma, jb, mb, jt.first.first, ma + mb)
                            * CgInt(jc, mc, jd, md, jt.first.first, ma + mb);
                        }
                        oaby6Matrix.coeffRef(alpha * orbitNumber + beta, gamma * orbitNumber + delta) += value;
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
        }
    }
    return index;
}


int kroneckerDelta(const int a, const int b) {
    return a == b ? 1 : 0;
}
