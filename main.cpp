#include <iostream>
#include <vector>

#include "CG.h"
#include "input.h"
#include "npa.h"

using namespace std;
int main() {

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

    vector<vector<pair<int, int>>> basesJZ;
    generateBasesJ(JisList, basesZ, basesJZ);

    vector<OrbitM*> orbitMProton, orbitMNeutron;
    initializeOrbitM(orbitProton, orbitMProton);
    initializeOrbitM(orbitNeutron, orbitMNeutron);

    vector<PairM*> pairZMs;
    vector<vector<int>> pairJMMap;

    initializePairM(pairZs, orbitProton, orbitMProton, pairZMs, pairJMMap);

    vector<vector<int>> basesZM0, basesZM1;

    generateMBases(pairZs, pairZMs, Z / 2, basesZ, pairJMMap, basesZM0, basesZM1);

    vector<vector<double>> transformMatrix;



    return 0;
}