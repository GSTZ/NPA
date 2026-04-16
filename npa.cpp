#include "CG.h"
#include "input.h"
#include "npa.h"

#include <string>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <set>

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


bool isSameVector(vector<int> a, vector<int> b) {
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    return a == b;
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
                    if (abs(v) > 1e-11) {
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


int overlapMScheme(const std::vector<PairM*>& pairMs,
                   const std::vector<std::vector<int>>& basesM,
                   const int orbitNumber,
                   Eigen::MatrixXd& overlapMap) {
    const int dim = static_cast<int>(basesM.size());
    overlapMap.resize(dim, dim);
    overlapMap.setZero();

    // rowStart[i] = 上三角中第 i 行的起始线性编号
    // 第 i 行包含 (i,i) 到 (i,dim-1)，共 dim-i 个元素
    std::vector<long long> rowStart(dim + 1, 0);
    for (int i = 0; i < dim; ++i) {
        rowStart[i + 1] = rowStart[i] + (dim - i);
    }

    const long long total = rowStart[dim];  // = dim*(dim+1)/2

    #pragma omp parallel for schedule(dynamic)
    for (long long k = 0; k < total; ++k) {
        // 找到满足 rowStart[i] <= k < rowStart[i+1] 的 i
        auto it = std::upper_bound(rowStart.begin(), rowStart.end(), k);
        int i = static_cast<int>(it - rowStart.begin()) - 1;

        long long offset = k - rowStart[i];
        int j = i + static_cast<int>(offset);

        const auto& basisMBra = basesM[i];
        const auto& basisMKet = basesM[j];

        double overlap = overlapMSchemeOne(pairMs, basisMBra, basisMKet, orbitNumber);

        overlapMap(i, j) = overlap;
        if (i != j) {
            overlapMap(j, i) = overlap;
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

    Eigen::SparseMatrix<double, Eigen::RowMajor> qbar;
    Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

    N1Bab(bra, ket, orbitNumber, qbar, Bab);

    Eigen::SparseMatrix<double, Eigen::RowMajor> pnBn = pairMs[basisMBra.back()]->pab * Bab;
    overlap = -2 * sparseTrace(pnBn);

    return overlap;
}


int N2Qaby6(const vector<PairNew*>& bra, const vector<PairNew*>& ket, const int orbitNumber,
    Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar) {

    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>> resultQ;

    generalMultiCommutator(bra, ket, resultQ);


    const int n = orbitNumber;
    const int qSize = static_cast<int>(resultQ.size());
    const int dim = n * n;
    const double one_sixth = 1.0 / 6.0;

    qbar.resize(dim, dim);

    // 24 个排列，以及对应奇偶性
    static const int perm[24][4] = {
        {0,1,2,3}, {0,1,3,2}, {0,2,1,3}, {0,2,3,1},
        {0,3,1,2}, {0,3,2,1}, {1,0,2,3}, {1,0,3,2},
        {1,2,0,3}, {1,2,3,0}, {1,3,0,2}, {1,3,2,0},
        {2,0,1,3}, {2,0,3,1}, {2,1,0,3}, {2,1,3,0},
        {2,3,0,1}, {2,3,1,0}, {3,0,1,2}, {3,0,2,1},
        {3,1,0,2}, {3,1,2,0}, {3,2,0,1}, {3,2,1,0}
    };

    static const int sign[24] = {
        +1,-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1,
        +1,-1,-1,+1,+1,-1,-1,+1,+1,-1,-1,+1
    };

    // 预缓存 lhs/rhs 指针，避免反复 resultQ[m].lhs()/rhs()
    std::vector<const Eigen::SparseMatrix<double>*> Ls(qSize);
    std::vector<const Eigen::SparseMatrix<double>*> Rs(qSize);
    for (int m = 0; m < qSize; ++m) {
        Ls[m] = std::vector<const Eigen::SparseMatrix<double> *>::value_type(&resultQ[m].lhs());
        Rs[m] = std::vector<const Eigen::SparseMatrix<double> *>::value_type(&resultQ[m].rhs());
    }

    // 每个线程一个 triplet 容器，避免锁竞争
    int numThreads = 1;
    #ifdef _OPENMP
    numThreads = omp_get_max_threads();
    #endif

    std::vector<std::vector<Eigen::Triplet<double>>> tripletsPerThread(numThreads);

    // 粗略预留一下空间：
    // 每个独立四元组最多产生 24 个非零
    // 独立四元组数量是 C(n,4)
    const long long numUnique =
        (n >= 4) ? (1LL * n * (n - 1) * (n - 2) * (n - 3) / 24) : 0LL;

    for (int t = 0; t < numThreads; ++t) {
        tripletsPerThread[t].reserve(static_cast<size_t>(numUnique * 24 / numThreads + 1024));
    }

    std::vector<std::array<int, 4>> combos;
    combos.reserve((n >= 4) ? (1LL * n * (n - 1) * (n - 2) * (n - 3) / 24) : 0);

    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            for (int c = b + 1; c < n; ++c) {
                for (int d = c + 1; d < n; ++d) {
                    combos.push_back({a, b, c, d});
                }
            }
        }
    }

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif

        auto& localTriplets = tripletsPerThread[tid];

        #pragma omp for schedule(dynamic)
        for (int idx = 0; idx < static_cast<int>(combos.size()); ++idx) {
            const int a = combos[idx][0];
            const int b = combos[idx][1];
            const int c = combos[idx][2];
            const int d = combos[idx][3];

            double value = 0.0;

            for (int m = 0; m < qSize; ++m) {
                const auto& L = *Ls[m];
                const auto& R = *Rs[m];

                const double qabcd = L.coeff(a, b) * R.coeff(c, d);
                const double qacbd = L.coeff(a, c) * R.coeff(b, d);
                const double qadbc = L.coeff(a, d) * R.coeff(b, c);
                const double qbcad = L.coeff(b, c) * R.coeff(a, d);
                const double qbdac = L.coeff(b, d) * R.coeff(a, c);
                const double qcdab = L.coeff(c, d) * R.coeff(a, b);

                value += (qabcd - qacbd + qadbc + qbcad - qbdac + qcdab) / 6.0;
            }

            if (std::abs(value) < 1e-11) {
                continue;
            }

            const int base[4] = {a, b, c, d};

            for (int p = 0; p < 24; ++p) {
                const int i = base[perm[p][0]];
                const int j = base[perm[p][1]];
                const int k = base[perm[p][2]];
                const int l = base[perm[p][3]];

                const int row = i * n + k;
                const int col = j * n + l;

                localTriplets.emplace_back(row, col, sign[p] * value);
            }
        }
    }

    // 合并 triplets
    std::vector<Eigen::Triplet<double>> allTriplets;
    size_t totalSize = 0;
    for (const auto& v : tripletsPerThread) {
        totalSize += v.size();
    }
    allTriplets.reserve(totalSize);

    for (auto& v : tripletsPerThread) {
        allTriplets.insert(allTriplets.end(), v.begin(), v.end());
    }

    // 批量构造稀疏矩阵；如果有重复项则自动求和
    qbar.setFromTriplets(allTriplets.begin(), allTriplets.end(),
                         [](const double& x, const double& y) { return x + y; });

    qbar.makeCompressed();

    return 0;
}


int N1Bab(const std::vector<PairNew*>& bra,
          const std::vector<PairNew*>& ket,
          const int orbitNumber,
          Eigen::SparseMatrix<double, Eigen::RowMajor>& qbar,
          Eigen::SparseMatrix<double, Eigen::RowMajor>& Bab)
{
    auto pN1 = bra.back()->pab;
    auto braNext = bra;
    braNext.pop_back();

    if (!(qbar.size() > 0)) N2Qaby6(braNext, ket, orbitNumber, qbar);

    const int n = orbitNumber;
    Bab.resize(n, n);
    Bab.setZero();

    // 先把 pN1 的非零元收集出来，避免后面反复扫
    struct NZ {
        int k;
        int l;
        double val;
    };
    std::vector<NZ> pList;
    pList.reserve(pN1.nonZeros());

    for (int k = 0; k < pN1.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(pN1, k); it; ++it) {
            pList.push_back({static_cast<int>(it.row()), static_cast<int>(it.col()), it.value()});
        }
    }

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(n) * n / 4 + 16);

    // 临时 dense 行缓存，减少对 Bab 的稀疏随机写
    std::vector<double> rowBuffer(n, 0.0);

    for (int i = 0; i < n; ++i) {
        std::fill(rowBuffer.begin(), rowBuffer.end(), 0.0);

        for (const auto& p : pList) {
            const int k = p.k;
            const int l = p.l;
            const double pv = p.val;

            const int qrow = i * n + k;

            for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(qbar, qrow); it; ++it) {
                const int col = it.col();

                // col = j*n + l2
                // 只保留 l2 == l 的项
                if (col % n == l) {
                    const int j = col / n;
                    rowBuffer[j] += 12.0 * pv * it.value();
                }
            }
        }

        for (int j = 0; j < n; ++j) {
            if (rowBuffer[j] != 0.0) {
                triplets.emplace_back(i, j, rowBuffer[j]);
            }
        }
    }

    Bab.setFromTriplets(triplets.begin(), triplets.end(),
                        [](double a, double b) { return a + b; });
    Bab.makeCompressed();

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

            if (norm_k_sq < 1e-8) {
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
        if (transformMatrix0.row(i).isZero(1e-11)) {
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
    Eigen::MatrixXd& gaby6Matrix, map<vector<int>, Eigen::SparseMatrix<double, Eigen::RowMajor>>& qbarMap) {
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
        Eigen::SparseMatrix<double, Eigen::RowMajor> qbar;
        Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

        auto bt2 = braTem;
        bt2.pop_back();
        sort(bt2.begin(), bt2.end());
        if (qbarMap.contains(bt2)) qbar = qbarMap[bt2];

        N1Bab(braNext, ketNext, orbitNumber, qbar, Bab);

        qbarMap[bt2] = qbar;

        for (int alpha = 0; alpha < orbitNumber; ++alpha) {
            for (int beta = 0; beta < orbitNumber; ++beta) {
                const auto pkab = 4.0 * Pk->pab.coeff(alpha, beta);
                gaby6Matrix.block(alpha * orbitNumber, beta * orbitNumber,
                                  orbitNumber, orbitNumber) += pkab * Bab;
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

            Eigen::SparseMatrix<double, Eigen::RowMajor> qbar;

            auto bt2 = braTem;
            sort(bt2.begin(), bt2.end());
            if (qbarMap.contains(bt2)) qbar = qbarMap[bt2];
            else {
                N2Qaby6(braNext, ketNext, orbitNumber, qbar);
                qbarMap[bt2] = qbar;
            }


            const int n = orbitNumber;
            const auto& PkMat = Pk->pab;
            const auto& PiMat = Pi->pab;

            for (int gamma = 0; gamma < n; ++gamma) {
                for (int delta = 0; delta < n; ++delta) {
                    std::vector<Eigen::Triplet<double>> triplets;
                    triplets.reserve(n * n / 10); // 按稀疏度估一个值

                    for (int alpha = 0; alpha < n; ++alpha) {
                        for (int beta = 0; beta < n; ++beta) {
                            double val = qbar.coeff(alpha * n + gamma, beta * n + delta);
                            if (std::abs(val) > 1e-11) {
                                triplets.emplace_back(alpha, beta, val);
                            }
                        }
                    }

                    if (triplets.empty()) continue;

                    Eigen::SparseMatrix<double, Eigen::RowMajor> qbarAB(n, n);
                    qbarAB.setFromTriplets(triplets.begin(), triplets.end());

                    Eigen::SparseMatrix<double, Eigen::RowMajor> Sab =
                        PkMat * qbarAB * PiMat + PiMat * qbarAB * PkMat;
                    Sab *= 96.0;

                    for (int k = 0; k < Sab.outerSize(); ++k) {
                        for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(Sab, k); it; ++it) {
                            gaby6Matrix(it.row() * n + gamma, it.col() * n + delta) += it.value();
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int fab(const vector<PairM*>& pairMs, const vector<int>& bra, const vector<int>& ket, const int orbitNumber,
    Eigen::MatrixXd& fabMatrix, map<vector<int>, Eigen::SparseMatrix<double, Eigen::RowMajor>>& qbarMap) {
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
        Eigen::SparseMatrix<double, Eigen::RowMajor> qbar;
        Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

        auto bt2 = braTem;
        bt2.pop_back();
        sort(bt2.begin(), bt2.end());
        if (qbarMap.contains(bt2)) qbar = qbarMap[bt2];

        N1Bab(braNext, ketNext, orbitNumber, qbar, Bab);
        fabMatrix += 4 * Pk->pab * Bab.transpose();

        qbarMap[bt2] = qbar;
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


int oaby6(const int orbitNumber,
          const std::vector<Orbit*>& orbitAC,
          const std::vector<OrbitM*>& orbitMAC,
          const std::vector<Orbit*>& orbitBD,
          const std::vector<OrbitM*>& orbitMBD,
          const std::map<TBMEJ, std::map<std::pair<int, int>, double>>& tbmeJMap,
          Eigen::SparseMatrix<double>& oaby6Matrix)
{
    const int n = orbitNumber;
    const int dim = n * n;

    oaby6Matrix.resize(dim, dim);
    oaby6Matrix.setZero();

    // 预提取 orbit 信息
    std::vector<int> aIndex(n), bIndex(n);
    std::vector<int> jA(n), jB(n);
    std::vector<int> mA(n), mB(n);

    for (int alpha = 0; alpha < n; ++alpha) {
        aIndex[alpha] = indexOrbitJM(orbitAC, orbitMAC[alpha]);
        jA[alpha] = orbitMAC[alpha]->j;
        mA[alpha] = orbitMAC[alpha]->m;
    }

    for (int beta = 0; beta < n; ++beta) {
        bIndex[beta] = indexOrbitJM(orbitBD, orbitMBD[beta]);
        jB[beta] = orbitMBD[beta]->j;
        mB[beta] = orbitMBD[beta]->m;
    }

    // 按 m 对 delta 分桶：给定 md，快速找到所有 delta
    std::unordered_map<int, std::vector<int>> mToDeltaList;
    mToDeltaList.reserve(n * 2);

    for (int delta = 0; delta < n; ++delta) {
        mToDeltaList[mB[delta]].push_back(delta);
    }

    int numThreads = 1;
    #ifdef _OPENMP
    numThreads = omp_get_max_threads();
    #endif

    std::vector<std::vector<Eigen::Triplet<double>>> tripletsPerThread(numThreads);

    #pragma omp parallel
    {
        int tid = 0;
        #ifdef _OPENMP
        tid = omp_get_thread_num();
        #endif

        auto& localTriplets = tripletsPerThread[tid];
        localTriplets.reserve(static_cast<size_t>(n) * n);

        #pragma omp for schedule(dynamic, 1)
        for (int alpha = 0; alpha < n; ++alpha) {
            const int a  = aIndex[alpha];
            const int ja = jA[alpha];
            const int ma = mA[alpha];

            for (int beta = 0; beta < n; ++beta) {
                const int b  = bIndex[beta];
                const int jb = jB[beta];
                const int mb = mB[beta];

                const int M = ma + mb;

                for (int gamma = 0; gamma < n; ++gamma) {
                    const int c  = aIndex[gamma];
                    const int jc = jA[gamma];
                    const int mc = mA[gamma];

                    const int mdTarget = M - mc;

                    auto foundDelta = mToDeltaList.find(mdTarget);
                    if (foundDelta == mToDeltaList.end()) {
                        continue;
                    }

                    const auto& deltaList = foundDelta->second;

                    for (int delta : deltaList) {
                        const int d  = bIndex[delta];
                        const int jd = jB[delta];
                        const int md = mB[delta];

                        // 这里按理已经满足 ma+mb == mc+md，但保守再检查一次也行
                        if (ma + mb != mc + md) {
                            continue;
                        }

                        const TBMEJ tj(a, b, c, d);
                        auto tbmeIt = tbmeJMap.find(tj);
                        if (tbmeIt == tbmeJMap.end()) {
                            continue;
                        }

                        const auto& JTList = tbmeIt->second;

                        const double normFactor =
                            std::sqrt((1.0 + kroneckerDelta(a, b)) *
                                      (1.0 + kroneckerDelta(c, d)));

                        double value = 0.0;

                        for (const auto& jt : JTList) {
                            const int J = jt.first.first;
                            const double tbmeVal = jt.second;

                            value += normFactor * tbmeVal
                                * CgInt(ja, ma, jb, mb, J, M)
                                * CgInt(jc, mc, jd, md, J, M);
                        }

                        if (std::abs(value) > 1e-11) {
                            const int row = alpha * n + gamma;
                            const int col = beta * n + delta;
                            localTriplets.emplace_back(row, col, value);
                        }
                    }
                }
            }
        }
    }

    std::vector<Eigen::Triplet<double>> allTriplets;
    size_t totalTriplets = 0;
    for (const auto& vec : tripletsPerThread) {
        totalTriplets += vec.size();
    }
    allTriplets.reserve(totalTriplets);

    for (auto& vec : tripletsPerThread) {
        allTriplets.insert(allTriplets.end(), vec.begin(), vec.end());
    }

    oaby6Matrix.setFromTriplets(allTriplets.begin(), allTriplets.end(),
                                [](const double& x, const double& y) { return x + y; });
    oaby6Matrix.makeCompressed();

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

    LambdaLanczos<double> engine(mv_mul, matrix.cols(), false, 5);
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


int calHamiltonianMatrix(const std::vector<PairM*>& pairMs,
    const std::vector<std::vector<int>>& basesM,
    const int orbitNumber,
    const std::vector<std::vector<Eigen::MatrixXd>>& fabMap,
    Eigen::SparseMatrix<double>& oaby6Matrix,
    Eigen::MatrixXd& qabMatrix,
    Eigen::MatrixXd& hamMatrix) {
    const int N = static_cast<int>(basesM.size());
    hamMatrix.resize(N, N);
    hamMatrix.setZero();

    // 预提取 oaby6Matrix 的非零元素
    struct SparseEntry {
        int row;
        int col;
        double value;
    };

    std::vector<SparseEntry> oaby6Entries;
    oaby6Entries.reserve(oaby6Matrix.nonZeros());

    for (int k = 0; k < oaby6Matrix.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(oaby6Matrix, k); it; ++it) {
            SparseEntry entry;
            entry.row = static_cast<int>(it.row());
            entry.col = static_cast<int>(it.col());
            entry.value = it.value();
            oaby6Entries.push_back(entry);
        }
    }


    #pragma omp parallel for collapse(1) schedule(dynamic)
    for (int j = 0; j < N; ++j) {

        const auto& basisMKet = basesM[j];
        map<vector<int>, Eigen::SparseMatrix<double, Eigen::RowMajor>> qbarMap;

        for (int i = j; i < N; ++i) {

            const auto& basisMBra = basesM[i];

            Eigen::MatrixXd gaby6Matrix;
            gaby6(pairMs, basisMBra, basisMKet, orbitNumber, gaby6Matrix, qbarMap);

            const Eigen::MatrixXd& fabMatrix = fabMap[i][j];

            double element = 0.0;

            // 稀疏 oaby6 部分
            for (const auto& entry : oaby6Entries) {
                element += entry.value * gaby6Matrix(entry.row, entry.col);
            }
            element /= 4.0;

            // fab 与 qab 的逐元素内积
            element += (fabMatrix.array() * qabMatrix.array()).sum();

            hamMatrix(i, j) = element;
            if (i != j) {
                hamMatrix(j, i) = element;
            }
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
                * (J + 1) * sixjSymbolsInt(ja, jc, J, jd, jb, t)
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
        //std::cerr << "Warning: Reducing rank from " << rank
        //          << " to max possible " << max_rank << std::endl;
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
    auto dim = ts.size();

    //#pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < rank; ++i) {
        for (int k = 0; k < dim; ++k) {
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
    #pragma omp parallel for collapse(2) schedule(dynamic)
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

    int dim1 = qabtMMatrixMap.size();
    int dim2 = qabtMMatrixMap[0].size();

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int m = 0; m < dim1; ++m) {
        for (int n = 0; n < dim2; ++n) {
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

    int dim1 = basesMBra.size();
    int dim2 = basesMKet.size();

    #pragma omp parallel for collapse(1) schedule(dynamic)
    for (int j = 0; j < dim2; j++) {
        const auto& basisMKet = basesMKet[j];
        map<vector<int>, Eigen::SparseMatrix<double, Eigen::RowMajor>> qbarMap;
        for (int i = 0; i < dim1; i++) {
            Eigen::MatrixXd fabMatrix;
            const auto& basisMBra = basesMBra[i];
            fab(pairMs, basisMBra, basisMKet, orbitNumber, fabMatrix, qbarMap);
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
           * sixjSymbolsInt(JNuPrime, JPiPrime, J, JPi, JNu, t)
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

    int dim1 = ts.size();

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < dim1; i++) {
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


// 阶乘函数 (支持分数)
double factorial1(double x) {
    if (x < 0) return 0; // 阶乘参数必须非负
    return std::tgamma(x + 1); // 使用伽玛函数计算阶乘
}

double factorial(int n) {
    if (n == 0) return 1;
    double result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}


// 求和部分的计算
double calculate_sum(int n, int n_prime, int l, int l_prime, int lambda) {
    double sum = 0.0;

    // 计算常量部分
    int l_l_prime_lambda = l + l_prime + lambda;

    // 遍历 q 的所有可能值
    for (int q = 0; q <= std::min(n, n_prime); ++q) {
        // 计算每一项中的分母因子
        int term1 = n - q;
        int term2 = n_prime - q;
        double term3 = q + (l_prime - l + lambda) / 2.0 - n;
        double term4 = q + (l - l_prime + lambda) / 2.0 - n_prime;

        // 确保所有参数非负，否则跳过该项
        if (term1 < 0 || term2 < 0 || term3 < 0 || term4 < 0) {
            continue;
        }

        // 计算分子部分 (l + l' + λ + 2q + 1)!!
        double numerator = double_factorial(l_l_prime_lambda + 2 * q + 1);

        // 计算分母部分
        double denominator = std::pow(2, q) * factorial1(q) *
                             factorial1(term1) * factorial1(term2) *
                             factorial1(term3) * factorial1(term4);

        // 累加到总和
        sum += numerator / denominator;
    }

    return sum;
}

// 主公式的计算
double calculate_integral(int n, int n_prime, int l, int l_prime, int lambda, double alpha) {
    // (-1)^(n - n')
    int sign_factor = (n - n_prime) % 2 == 0 ? 1 : -1;

    // 开平方项
    double sqrt_factor = std::sqrt(
        (factorial1(n) * factorial1(n_prime) * std::pow(2, n + n_prime - lambda)) /
        (double_factorial(2 * n + 2 * l + 1) * double_factorial(2 * n_prime + 2 * l_prime + 1))
    );

    // 对称性项 (1 + (-1)^(l + l' + λ)) / 2
    int symmetry_factor = (l + l_prime + lambda) % 2 == 0 ? 1 : 0;

    // 阶乘部分
    double factorial_part = factorial1((l_prime - l + lambda) / 2.0) *
                            factorial1((l - l_prime + lambda) / 2.0);

    // 求和部分
    double summation = calculate_sum(n, n_prime, l, l_prime, lambda);

    // 最终结果
    return sign_factor * sqrt_factor * symmetry_factor * factorial_part * summation / std::pow(alpha, lambda);
}

double computeQ(int j11, int j22, double tt, int n11, int l11, int n22, int l22,double alpha) {
    double j1=j11/2.0;
    double j2=j22/2.0;
    double t=tt/2.0;
    int n1=n11/2;
    int n2=n22/2;
    int l1=l11/2;
    int l2=l22/2;
    // 计算前面的符号部分
    double prefactor = std::pow(-1, j1 + t - 0.5);

    // 计算根号部分
    double sqrtPart = std::sqrt((2 * j1 + 1) * (2 * j2 + 1) / (4 * M_PI * (2 * t + 1)));

    // 计算 Clebsch-Gordan 系数
    double clebsch = CgInt(j11, j22, tt, 1, -1,0);
    // double alphaa=4.33158;

    // 计算径向积分
    double radialIntegral = calculate_integral(n1, n2, l1, l2, t, alpha);

    // 计算最终结果
    double result = prefactor * sqrtPart * clebsch * (1 + std::pow(-1, l1 + l2 + t)) / 2.0 * radialIntegral;
    if (t==0)
    {
        if (j11==j22)
        {
            result=sqrt(2.0*j1+1)/2;
        }
        else
        {
            return 0;
        }
    }
    return result;
}


// ============================================================
//  广义 Laguerre 多项式 L_n^alpha(x)
//  使用三项递推关系（数值稳定）
//    L_0^a(x) = 1
//    L_1^a(x) = 1 + a - x
//    (n+1)*L_{n+1} = (2n+1+a-x)*L_n - (n+a)*L_{n-1}
// ============================================================
double laguerre(int n, double alpha, double x) {
    if (n == 0) return 1.0;
    double L0 = 1.0;
    double L1 = 1.0 + alpha - x;
    if (n == 1) return L1;
    double Lk = 0.0;
    for (int k = 1; k < n; ++k) {
        Lk = ((2*k + 1 + alpha - x) * L1 - (k + alpha) * L0) / (k + 1);
        L0 = L1;
        L1 = Lk;
    }
    return Lk;
}

// ============================================================
//  谐振子径向波函数 R_{nl}(r)（b=1 无量纲）
//  R_{nl}(r) = N_{nl} * r^l * exp(-r^2/2) * L_n^{l+1/2}(r^2)
//  满足归一化：∫ |R_{nl}|^2 r^2 dr = 1
// ============================================================
double R_nl(int n, int l, double r) {
    // 归一化系数 N = sqrt(2*n! / Gamma(n+l+3/2))
    double norm = std::sqrt(2.0 * std::tgamma(n + 1) /
                            std::tgamma(n + l + 1.5));
    double x    = r * r;
    double lag  = laguerre(n, l + 0.5, x);
    return norm * std::pow(r, l) * std::exp(-x / 2.0) * lag;
}

// ============================================================
//  Gauss-Legendre 积分节点和权重（N 点）
//  区间 [a, b] 上的积分
//  使用迭代法求 Legendre 多项式零点
// ============================================================
void gauss_legendre(int N, double a, double b,
                    std::vector<double>& x,
                    std::vector<double>& w) {
    x.resize(N);
    w.resize(N);
    const double PI = std::acos(-1.0);

    for (int i = 0; i < (N + 1) / 2; ++i) {
        // 初始猜测
        double xi = std::cos(PI * (i + 0.75) / (N + 0.5));
        double dP, P;
        // Newton 迭代
        for (int iter = 0; iter < 100; ++iter) {
            double P0 = 1.0, P1 = xi;
            for (int j = 2; j <= N; ++j) {
                P  = ((2*j - 1) * xi * P1 - (j - 1) * P0) / j;
                P0 = P1;
                P1 = P;
            }
            P  = P1;
            dP = N * (P0 - xi * P1) / (1.0 - xi * xi);
            double dx = P / dP;
            xi -= dx;
            if (std::abs(dx) < 1e-15) break;
        }
        // 计算权重
        double P0 = 1.0, P1 = xi;
        for (int j = 2; j <= N; ++j) {
            P  = ((2*j - 1) * xi * P1 - (j - 1) * P0) / j;
            P0 = P1;
            P1 = P;
        }
        dP = N * (P0 - xi * P1) / (1.0 - xi * xi);

        double wi = 2.0 / ((1.0 - xi*xi) * dP * dP);

        // 映射到 [a, b]
        x[i]       = 0.5 * (b - a) * (-xi) + 0.5 * (b + a);
        x[N-1-i]   = 0.5 * (b - a) * ( xi) + 0.5 * (b + a);
        w[i]       = 0.5 * (b - a) * wi;
        w[N-1-i]   = w[i];
    }
}

// ============================================================
//  主积分函数：<n1,l1 | r^lambda | n2,l2>
//
//  策略：
//    - 积分区间 [0, r_max]，r_max 自适应取 max_r(n,l)
//    - 使用 Gauss-Legendre 积分（N_pts 点）
//    - 被积函数：R_{n1,l1}(r) * r^lambda * R_{n2,l2}(r) * r^2
// ============================================================
double radial_matrix_element_numerical(
    int n1, int l1, int n2, int l2, int lam,
    double b,           // 谐振子长度参数（fm）
    int    N_pts = 200  // 积分点数，默认 200 点
) {
    // 宇称选择定则
    if ((l1 + l2 + lam) % 2 != 0) return 0.0;

    // 自适应确定积分上限（波函数在 r_max 处已充分衰减）
    // 经验公式：r_max ≈ (2*max_n + max_l + 6) * sqrt(2) * b
    int    n_max  = std::max(n1, n2);
    int    l_max  = std::max(l1, l2);
    double r_max  = (2.0 * n_max + l_max + 8.0) * std::sqrt(2.0);
    // 注意：无量纲化后 b=1，最后乘 b^lambda

    // 获取 Gauss-Legendre 节点和权重
    std::vector<double> r_pts, wts;
    gauss_legendre(N_pts, 0.0, r_max, r_pts, wts);

    // 数值积分
    double integral = 0.0;
    for (int i = 0; i < N_pts; ++i) {
        double r  = r_pts[i];
        double f  = R_nl(n1, l1, r)
                  * std::pow(r, lam)
                  * R_nl(n2, l2, r)
                  * r * r;           // r^2 来自球坐标体积元
        integral += wts[i] * f;
    }

    // 乘以 b^lambda（量纲还原）
    return integral * std::pow(b, static_cast<double>(lam)) / std::pow(b, static_cast<double>(2 * lam));
}


double calculateBE2(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& braBasesJPN,
    const vector<vector<int>>& ketBasesJPN,
    const int t,
    const State& initial,
    const State& final,
    const double ePi,
    const double eNu,
    const Eigen::MatrixXd& reducedQopJP,
    const Eigen::MatrixXd& reducedQopJN) {
    double result = 0.0;
    double tFi = calTfi(basesJP, basesJN, braBasesJPN, ketBasesJPN, t, initial, final, ePi, eNu, reducedQopJP, reducedQopJN);
    result = (final.J + 1.0) / (initial.J + 1.0) * tFi * tFi;
    return result;
}


int generateQj1j2J(const vector<Orbit*>& orbits,
    const int t,
    const double alpha,
    Eigen::MatrixXd& qj1j2J) {
    auto orbitNumber = orbits.size();
    qj1j2J.resize(orbitNumber, orbitNumber);
    qj1j2J.setZero();

    for (int i = 0; i < orbitNumber; ++i) {
        for (int j = 0; j < orbitNumber; ++j) {
            const int n1 = orbits[i]->n;
            const int l1 = orbits[i]->l;
            const int j1 = orbits[i]->j;   // 2*j

            const int n2 = orbits[j]->n;
            const int l2 = orbits[j]->l;
            const int j2 = orbits[j]->j;   // 2*j

            const double phase = ((((j1 - 1) / 2) % 2) == 0) ? 1.0 : -1.0;
            const double phase2 = (1.0 + std::pow(-1, l1 + l2)) / 2.0;
            const double jhat1 = std::sqrt(static_cast<double>(j1 + 1));
            const double jhat2 = std::sqrt(static_cast<double>(j2 + 1));
            const double cg = CgInt(j1, 1, j2, -1, t, 0);

            const double radial1 = calculate_integral(n1, n2, l1, l2, t / 2, alpha);
            const double radial2 = radial_matrix_element_numerical(n1, l1, n2, l2, t / 2, alpha, 200);

            if (abs(radial1 - radial2) > 1e-6) { cout << radial1 << "  " << radial2 << endl; }

            qj1j2J(i, j) =
                phase * phase2 * jhat1 * jhat2 / std::sqrt(20.0 * M_PI) * cg * radial2;
        }
    }
    return 0;
}



int qj1j2JtoM(const vector<OrbitM*>& orbitMs,
    const vector<Orbit*>& orbits,
    const Eigen::MatrixXd& qj1j2J,
    const int t,
    Eigen::MatrixXd& qj1j2M) {
    auto orbitNumber = orbitMs.size();
    qj1j2M.resize(orbitNumber, orbitNumber);
    qj1j2M.setZero();
    for (int i = 0; i < orbitNumber; ++i) {
        auto bra = orbitMs[i];
        auto JAlpha = indexOrbitJM(orbits, bra);
        for (int j = 0; j < orbitNumber; ++j) {
            auto ket = orbitMs[j];
            auto JBeta = indexOrbitJM(orbits, ket);
            qj1j2M(i, j) = std::pow(-1, (ket->j + ket->m) / 2.0)
            * CgInt(bra->j, bra->m, ket->j, -ket->m, t, bra->m - ket->m)
            * qj1j2J(JAlpha, JBeta);
        }
    }
    return 0;
}


int generateQOperatorM(const vector<vector<int>>& braBasesM,
    const vector<vector<int>>& ketBasesM,
    const std::vector<std::vector<Eigen::MatrixXd>>& fabMap,
    const Eigen::MatrixXd& qj1j2M,
    Eigen::MatrixXd& Qop) {
    const int dimBra = static_cast<int>(braBasesM.size());
    const int dimKet = static_cast<int>(ketBasesM.size());

    Qop.resize(dimBra, dimKet);
    Qop.setZero();

    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < dimBra; ++i) {
        for (int j = 0; j < dimKet; ++j) {
            Qop(i, j) = (fabMap[i][j].cwiseProduct(qj1j2M)).sum();
        }
    }
    return 0;
}


int generateReducedOrthogonalQopJ(
    const vector<vector<pair<int, int>>>& basesJ,
    const Eigen::MatrixXd& QopJ00,
    const Eigen::MatrixXd& QopJ01,
    const Eigen::MatrixXd& orthogonalBasis,
    const int t,
    Eigen::MatrixXd& reducedJMatrix) {
    int dim1 = basesJ.size();

    reducedJMatrix.resize(dim1, dim1);
    reducedJMatrix.setZero();
    for (int k = 0; k < dim1; k++) {
        for (int l = 0; l < dim1; l++) {
            int J = basesJ[k].back().second;
            int JPrime = basesJ[l].back().second;
            if (isTriangle(JPrime, J, t)) {
                if ((J + JPrime + t) % 4 == 0) {
                    reducedJMatrix(k, l) = QopJ00(k, l) / CgInt(JPrime, 0, t, 0, J, 0);
                } else {
                    reducedJMatrix(k, l) = QopJ01(k, l) / CgInt(JPrime, 2, t, -2, J, 0);
                }
            }
        }
    }
    reducedJMatrix = orthogonalBasis * reducedJMatrix * orthogonalBasis.transpose();

    return 0;
}


double calTfi(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& braBasesJPN,
    const vector<vector<int>>& ketBasesJPN,
    const int t,
    const State& initial,
    const State& final,
    const double ePi,
    const double eNu,
    const Eigen::MatrixXd& ReducedQopJP,
    const Eigen::MatrixXd& ReducedQopJN) {

    double result = 0.0;

    int dim1 = final.waveFunction.size();
    int dim2 = initial.waveFunction.size();

    for (int i1 = 0; i1 < dim1; ++i1) {
        for (int i2 = 0; i2 < dim2; ++i2) {
            const auto& bra = braBasesJPN[i1];
            const auto& ket = ketBasesJPN[i2];
            double finalFactor = final.waveFunction[i1];
            double initialFactor = initial.waveFunction[i2];

            result += finalFactor * (ePi * JfQtJi(bra, ket, basesJP, basesJN, t, 0, ReducedQopJP, ReducedQopJN)
                + eNu * JfQtJi(bra, ket, basesJP, basesJN, t, 1, ReducedQopJP, ReducedQopJN)) * initialFactor;
        }
    }

    return result;
}

double JfQtJi(const vector<int>& bra, const vector<int>& ket, // 0 Pi, 1 Nu, 2 J, 3 M
    const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const int t,
    const int PiNu, // 0 Pi, 1 Nu
    const Eigen::MatrixXd& ReducedQopJP,
    const Eigen::MatrixXd& ReducedQopJN) {
    double result = 0.0;
    int JPiPrime = basesJP[bra[0]].back().second;
    int JNuPrime = basesJN[bra[1]].back().second;
    int JPi = basesJP[ket[0]].back().second;
    int JNu = basesJN[ket[1]].back().second;
    int JFinal = bra[2];
    int JInitial = ket[2];

    int braJPIndex = bra[0];
    int ketJPIndex = ket[0];
    int braJNIndex = bra[1];
    int ketJNIndex = ket[1];

    if (PiNu == 0) {
        if (braJNIndex == ketJNIndex) result = std::pow(-1, (JNuPrime + JPi + JFinal + t) / 2.0)
            * sqrt(JInitial + 1.0) * sqrt(JPiPrime + 1.0)
            * sixjSymbolsInt(JNuPrime, JPi, JInitial, t, JFinal, JPiPrime)
            * ReducedQopJP(braJPIndex, ketJPIndex);
    } else {
        if (braJPIndex == ketJPIndex) result = std::pow(-1, (JNu - JNuPrime + JFinal - JInitial + JPiPrime + JNu + JFinal + t) / 2.0)
            * sqrt(JInitial + 1.0) * sqrt(JNuPrime + 1.0)
            * sixjSymbolsInt(JPiPrime, JNu, JInitial, t, JFinal, JNuPrime)
            * ReducedQopJN(braJNIndex, ketJNIndex);
    }
    return result;
}


int generateMj1j2J(const vector<Orbit*>& orbits,
    const int t,
    const double alpha,
    const double gs,
    const double gl,
    const double mun,
    Eigen::MatrixXd& mj1j2J) {
    auto orbitNumber = orbits.size();
    mj1j2J.resize(orbitNumber, orbitNumber);
    mj1j2J.setZero();

    for (int i = 0; i < orbitNumber; ++i) {
        for (int j = 0; j < orbitNumber; ++j) {
            const int n1 = orbits[i]->n;
            const int l1 = orbits[i]->l;
            const int j1 = orbits[i]->j;   // 2*j

            const int n2 = orbits[j]->n;
            const int l2 = orbits[j]->l;
            const int j2 = orbits[j]->j;   // 2*j

            if (l1 == l2 and n1 == n2) {

                // const double jhat1 = std::sqrt(static_cast<double>(j1 + 1));
                // const double jhat2 = std::sqrt(static_cast<double>(j2 + 1));
                //
                // double term1 = gl * std::pow(-1.0, l1 + 0.5 + j2 / 2.0)
                //             * sqrt(l1 * (l1 + 1.0) / 3.0)
                //             * jhat1 * jhat2 * sqrt(2 * l1 + 1.0)
                //             * sixjSymbolsInt(j1, j2, 2, 2 * l1, 2 * l1, 1);
                //
                // double term2 = gs * std::pow(-1.0, l1 + 0.5 + j1 / 2.0)
                //             * sqrt(0.5) * jhat1 * jhat2
                //             * sixjSymbolsInt(j1, j2, 2, 1, 1, 2 * l1);
                //
                // mj1j2J(i, j) = 1 * (term1 + term2);


                const double jhat1 = std::sqrt(static_cast<double>(j1 + 1));
                const double jhat2 = std::sqrt(static_cast<double>(j2 + 1));

                 const double factor1 = mun * sqrt(t / 2.0 * (t + 1.0));
                 const double term1 = 2.0 / (t / 2.0 + 1.0) * gl
                             * std::pow(-1.0, j1 / 2.0 + j2 / 2.0 + l1 + l2 + t / 2.0 + 1.0)
                             * (1.0 + std::pow(-1.0, l1 + l2 + t / 2.0 + 1.0)) / 2.0
                             * sqrt(j2 / 2.0 * (j2 / 2.0 + 1) / (4 * M_PI))
                             * jhat1 * jhat2 * sqrt(t - 1.0)
                             * CgInt(j1, 1, t - 2, 0, j2, 1)
                             * sixjSymbolsInt(j1, t, j2, 2, j2, t - 2);
                 const double term2 = (gs - 2.0 / (t / 2.0 + 1.0) * gl)
                             * sqrt(1.0 / (8.0 * M_PI * t))
                             * std::pow(-1.0, t / 2.0)
                             * (1.0 + std::pow(-1.0, l1 + l2 + t / 2.0 + 1.0)) / 2.0
                             * jhat1 / sqrt(t + 1.0)
                             * (t / 2.0 - (std::pow(-1.0, (j1 + 1.0) / 2.0 + l1) * (j1 + 1.0) + std::pow(-1.0, (j2 + 1.0) / 2.0 + l2) * (j2 + 1.0)) / 2.0)
                             * CgInt(j1, 1, t, 0, j2, 1);

                const double radial1 = calculate_integral(n1, n2, l1, l2, t / 2 - 1, alpha);
                //const double radial2 = radial_matrix_element_numerical(n1, l1, n2, l2, t / 2 - 1, alpha, 200);

                mj1j2J(i, j) = factor1 * (term1 + term2) * radial1;
            }

            //const double jhat1 = std::sqrt(static_cast<double>(j1 + 1));
            //const double jhat2 = std::sqrt(static_cast<double>(j2 + 1));

            // const double factor1 = mun * sqrt(t / 2.0 * (t + 1.0));
            // const double term1 = 2.0 / (t / 2.0 + 1.0) * gl
            //             * std::pow(-1.0, j1 / 2.0 + j2 / 2.0 + l1 + l2 + t / 2.0 + 1.0)
            //             * (1.0 + std::pow(-1.0, l1 + l2 + t / 2.0 + 1.0)) / 2.0
            //             * sqrt(j2 / 2.0 * (j2 / 2.0 + 1) / (4 * M_PI))
            //             * jhat1 * jhat2 * sqrt(t - 1.0)
            //             * CgInt(j1, 1, t - 2, 0, j2, 1)
            //             * sixjSymbolsInt(j1, t, j2, 2, j2, t - 2);
            // const double term2 = (gs - 2.0 / (t / 2.0 + 1.0) * gl)
            //             * sqrt(1.0 / (8.0 * M_PI * t))
            //             * std::pow(-1.0, t / 2.0)
            //             * (1.0 + std::pow(-1.0, l1 + l2 + t / 2.0 + 1.0)) / 2.0
            //             * jhat1 / sqrt(t + 1.0)
            //             * (t / 2.0 - (std::pow(-1.0, (j1 + 1.0) / 2.0 + l1) * (j1 + 1.0) + std::pow(-1.0, (j2 + 1.0) / 2.0 + l2) * (j2 + 1.0)) / 2.0)
            //             * CgInt(j1, 1, t, 0, j2, 1);

            //const double radial1 = calculate_integral(n1, n2, l1, l2, t / 2 - 1, alpha);
            //const double radial2 = radial_matrix_element_numerical(n1, l1, n2, l2, t / 2 - 1, alpha, 200);

            //mj1j2J(i, j) = factor1 * (term1 + term2) * radial2;


        }
    }
    return 0;
}


double calMfi(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& braBasesJPN,
    const vector<vector<int>>& ketBasesJPN,
    const int t,
    const State& initial,
    const State& final,
    const Eigen::MatrixXd& ReducedMopJP,
    const Eigen::MatrixXd& ReducedMopJN) {

    double result = 0.0;

    int dim1 = final.waveFunction.size();
    int dim2 = initial.waveFunction.size();

    for (int i1 = 0; i1 < dim1; ++i1) {
        for (int i2 = 0; i2 < dim2; ++i2) {
            const auto& bra = braBasesJPN[i1];
            const auto& ket = ketBasesJPN[i2];
            double finalFactor = final.waveFunction[i1];
            double initialFactor = initial.waveFunction[i2];

            result += finalFactor * (JfQtJi(bra, ket, basesJP, basesJN, t, 0, ReducedMopJP, ReducedMopJN)
                + JfQtJi(bra, ket, basesJP, basesJN, t, 1, ReducedMopJP, ReducedMopJN)) * initialFactor;
        }
    }

    return result;
}


double calculateBM1(const vector<vector<pair<int, int>>>& basesJP,
    const vector<vector<pair<int, int>>>& basesJN,
    const vector<vector<int>>& braBasesJPN,
    const vector<vector<int>>& ketBasesJPN,
    const int t,
    const State& initial,
    const State& final,
    const Eigen::MatrixXd& reducedQopJP,
    const Eigen::MatrixXd& reducedQopJN) {
    double result = 0.0;
    double mFi = calMfi(basesJP, basesJN, braBasesJPN, ketBasesJPN, t, initial, final, reducedQopJP, reducedQopJN);
    result = (final.J + 1.0) / (initial.J + 1.0) * mFi * mFi;
    return result;
}