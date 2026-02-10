#include "CG.h"
#include "input.h"
#include "npa.h"

#include <string>
#include <vector>

using namespace std;

template <typename T, int Options>
T sparseTrace(const Eigen::SparseMatrix<T, Options>& mat) {
    T trace = T(0);
    for (Eigen::Index i = 0; i < std::min(mat.rows(), mat.cols()); ++i) {
        trace += mat.coeff(i, i);
    }
    return trace;
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
    const vector<vector<int>>& bases, vector<vector<int>>& basesM, const vector<vector<int>>& pairJMMap) {
    for (const auto& basis : bases) {
        vector<int> basisM;
        generateMBasis(pairs, pairMs, pairNumber, basis, basisM, pairJMMap, basesM);
    }
    return 0;
}

int generateMBasis(const vector<Pair*>& pairs, const vector<PairM*>& pairMs, const int& pairNumber,
    const vector<int>& basis, const vector<int>& basisM, const vector<vector<int>>& pairJMMap,
    vector<vector<int>>& basesM) {
    int index = basisM.size();
    if (index == pairNumber) {
        if (judgeMBasis(pairMs, basisM)) basesM.push_back(basisM);
    } else {
        auto pair = pairs[basis[index]];
        for (const auto pmIndex : pairJMMap[pair->index]) {
            auto pairM = pairMs[pmIndex];
            vector<int> basisMNext = basisM;
            basisMNext.push_back(pairM->index);
            generateMBasis(pairs, pairMs, pairNumber, basis, basisMNext, pairJMMap, basesM);
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


bool judgeMBasis(const vector<PairM*>& pairMs, const vector<int>& basisM) {
    bool flag = false;
    int M = 0;
    for (const auto& bmi : basisM) {
        M += pairMs[bmi]->m;
    }
    if (M == 0 or M == 2) flag = true;
    return flag;
}


int overlapMScheme();

int overlapMSchemeOne(const vector<PairM*>& pairMs, const vector<int>& basisMBra,
    const vector<int>& basisMKet, const int orbitNumber) {
    double overlap = 0.0;
    vector<Eigen::Product<Eigen::SparseMatrix<double, Eigen::RowMajor>, Eigen::SparseMatrix<double, Eigen::RowMajor>,
    Eigen::AliasFreeProduct>> resultQ;
    vector<PairNew*> bra, ket;
    // bra and ket have the same length.
    for (int i = 0; i < basisMBra.size() - 2; i++) {
        auto pNewBra = new PairNew();
        pNewBra->pab = pairMs[basisMBra[i]]->pab;
        bra.push_back(pNewBra);
    }

    for (int i = 0; i < basisMKet.size(); i++) {
        auto pNewKet = new PairNew();
        pNewKet->pab = pairMs[basisMKet[i]]->pab;
        ket.push_back(pNewKet);
    }
    generalMultiCommutator(bra, ket, resultQ);

    //vector<vector<Eigen::SparseMatrix<double, Eigen::RowMajor>>> qbar;
    Eigen::SparseMatrix<double, Eigen::RowMajor> qbar(orbitNumber ^ 2, orbitNumber ^ 2);


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
                            + resultQ[m].lhs().coeff(j, l) * resultQ[m].rhs().coeff(i, k)) / 6;
                    }
                }
            }
        }
    }

    Eigen::SparseMatrix<double, Eigen::RowMajor> Bab(orbitNumber, orbitNumber);

    for (int i = 0; i < orbitNumber; i++) {
        for (int j = 0; j < orbitNumber; j++) {
            for (int k = 0; k < orbitNumber; k++) {
                for (int l = 0; l < orbitNumber; l++) {
                    // to do
                }
            }
        }
    }

    Bab *= 12;

    for (int i = 0; i < basisMBra.size(); i++) {
        delete bra[i];
        delete ket[i];
    }
    return overlap;
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
