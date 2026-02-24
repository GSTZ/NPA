#include <iostream>
#include <vector>

#include "CG.h"
#include "input.h"
#include "npa.h"

using namespace std;
int main() {
    // int j1, j2, j3, m1, m2, m3;
    //
    // for (j1 = 0; j1 < 6; j1++) {
    //     for (j2 = 0; j2 < 6; j2++) {
    //         for (j3 = 0; j3 < 6; j3++) {
    //             for (m1 = 0; m1 < 6; m1++) {
    //                 for (m2 = 0; m2 < 6; m2++) {
    //                     for (m3 = 0; m3 < 6; m3++) {
    //                         double sixj = sixjSymbols(static_cast<double>(j1) / 2, static_cast<double>(j2) / 2,
    //                             static_cast<double>(j3) / 2, static_cast<double>(m1) / 2,
    //                             static_cast<double>(m2) / 2, static_cast<double>(m3) / 2);
    //                         if (abs(sixj) > 0.000000001) {
    //                             cout << sixj << endl;
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    //double sixj = sixjSymbols(0.5, 1, 1.5, 2, 2.5, 2);
    //cout << sixj << endl;
    string spsFile, pairFile;
    int Z, N, interactionNumber, orbitProtonNumber, orbitNeutronNumber, pairTypes, limitNumber;
    vector<Orbit*> orbitProton, orbitNeutron;
    vector<InteractionFile*> interactionFiles;
    vector<Pair*> pairZs;
    vector<PairLimit*> limits;
    readInput(Z, N, spsFile, pairFile, interactionNumber, interactionFiles);
    readSps(spsFile, orbitProtonNumber, orbitNeutronNumber, orbitProton, orbitNeutron);
    readPair(pairFile, pairTypes, pairZs, limitNumber, limits);

    for (int i = 0; i < pairZs.size(); i++) {
        pairZs[i]->yabr.setZero(orbitProtonNumber, orbitNeutronNumber);
    }

    vector<int> basis;
    vector<vector<int>> basesZ;

    generateBasis(pairZs, Z / 2, basis, basesZ, limits);

    vector<vector<vector<int>>> JisList;
    generateJi(pairZs, basesZ, JisList);

    vector<OrbitM*> orbitMProton, orbitMNeutron;
    initializeOrbitM(orbitProton, orbitMProton);
    initializeOrbitM(orbitNeutron, orbitMNeutron);

    vector<PairM*> pairZMs;
    vector<vector<int>> pairJMMap;

    initializePairM(pairZs, orbitProton, orbitMProton, pairZMs, pairJMMap);

    vector<vector<int>> basesZM;

    generateMBases(pairZs, pairZMs, Z / 2, basesZ, pairJMMap, basesZM);


    return 0;
}