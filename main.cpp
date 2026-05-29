#include <iostream>
#include <vector>

#include "CG.h"
#include "input.h"
#include "npa.h"
#include <chrono>

using namespace std;
int main() {

    auto start = std::chrono::high_resolution_clock::now();

    // single particle orbit file, collective pair file, structure coefficient file
    string spsFile, pairFile, scFileP, scFileN, interactionType;
    int Z, N, orbitPNumber, orbitNNumber, pairTypes, limitNumber;
    vector<Orbit*> orbitP, orbitN, nonsense;
    vector<vector<InteractionFile*>> interactionFiles;
    vector<Pair*> pairPs, pairNs;
    vector<PairLimit*> limits;
    vector<int> Js;

    // a, b, c, d, J, T, value
    map<TBMEJ, map<pair<int, int>, double>> tbmeJMap, tbmePN, tbmePP, tbmeNN;

    // read input files
    readInput(Z, N, spsFile, pairFile, scFileP, scFileN, interactionType, interactionFiles, Js);
    readSps(spsFile, orbitPNumber, orbitNNumber, orbitP, orbitN);

    readPair(pairFile, pairTypes, pairPs, limitNumber, limits);
    readPair(pairFile, pairTypes, pairNs, limitNumber, limits);

    readStructureCoefficient(scFileP, pairPs, pairTypes, orbitPNumber);
    readStructureCoefficient(scFileN, pairNs, pairTypes, orbitNNumber);

    if (interactionType == "pn") {
        for (auto& ifs : interactionFiles[0]) {
            readInteraction(ifs->filename, orbitP, tbmePP);
            interactionAntiSymmetric(tbmePP, orbitP);
        }

        for (auto& ifs : interactionFiles[1]) {
            readInteraction(ifs->filename, orbitN, tbmeNN);
            interactionAntiSymmetric(tbmeNN, orbitN);
        }

        for (auto& ifs : interactionFiles[2]) {
            readInteraction(ifs->filename, nonsense, tbmePN);
            interactionAntiSymmetricPN(tbmePN);
        }
    }

    //readInteraction(interactionFiles[0]->filename, orbitP, orbitN, tbmeJMap);
    //readInteraction(interactionFiles[1]->filename, orbitP, orbitN, tbmePN);

    bool ZOdd = Z % 2 != 0;
    bool NOdd = N % 2 != 0;

    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "read file duration: " << duration.count() << " ms" << std::endl;

    start = end;


    vector<int> basis;
    vector<vector<int>> basesP, basesN;
    vector<vector<vector<int>>> JisListP, JisListN;

    // In the order of J, pair<collective pair, Ji>.
    vector<vector<pair<int, int>>> basesJP, basesJN;

    // start and end position of one J.
    vector<pair<int, int>> blockJP, blockJN;

    vector<OrbitM*> orbitMProton, orbitMNeutron;

    vector<PairM*> pairPMs, pairNMs;
    vector<vector<int>> pairJMMapP, pairJMMapN;
    vector<vector<int>> basesPM0, basesPM1, basesNM0, basesNM1;

    // generate orbits and bases
    vector<A0*> A0sP, A0sN;
    if (ZOdd) generateA0s(orbitP, A0sP);
    if (NOdd) generateA0s(orbitN, A0sN);

    vector<P0*> P0sP, P0sN;
    if (ZOdd) generateP0s(A0sP, P0sP);
    if (NOdd) generateP0s(A0sN, P0sN);

    generateBasis(pairPs, Z / 2, orbitP.size(), ZOdd, basis, basesP, limits);
    generateBasis(pairNs, N / 2, orbitN.size(), NOdd, basis, basesN, limits);

    generateJi(pairPs, A0sP, basesP, JisListP);
    generateJi(pairNs, A0sN, basesN, JisListN);

    generateBasesJ(JisListP, basesP, basesJP);
    generateBasesJ(JisListN, basesN, basesJN);

    //buildBlockJ(basesJP, blockJP);
    //buildBlockJ(basesJN, blockJN);

    initializeOrbitM(orbitP, orbitMProton);
    initializeOrbitM(orbitN, orbitMNeutron);

    initializePairM(pairPs, orbitP, orbitMProton, pairPMs, pairJMMapP);
    initializePairM(pairNs, orbitN, orbitMNeutron, pairNMs, pairJMMapN);


    // prepare new pairs
    vector<vector<PairNew*>> preCalPairPMs, preCalPairNMs;
    const int levelEndP = 1;
    const int levelEndN = 1;

    preCalculateNewPairs(0, levelEndP, pairPMs, preCalPairPMs);
    preCalculateNewPairs(0, levelEndN, pairNMs, preCalPairNMs);


    generateMBases(pairPs, pairPMs, Z / 2, A0sP, P0sP, basesP, pairJMMapP, basesPM0, basesPM1);
    generateMBases(pairNs, pairNMs, N / 2, A0sN, P0sN, basesN, pairJMMapN, basesNM0, basesNM1);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "generate bases duration: " << duration.count() << " ms" << std::endl;

    start = end;


    Eigen::MatrixXd overlapMP00, overlapJP, orthogonalBasisP, overlapMN00, overlapJN, orthogonalBasisN, overlapMP11,
    overlapMN11;

    // calculate M-scheme overlap
    overlapMScheme(pairPMs, basesPM0, ZOdd, P0sP, orbitMProton.size(), levelEndP, preCalPairPMs, overlapMP00);
    overlapMScheme(pairNMs, basesNM0, NOdd, P0sN, orbitMNeutron.size(), levelEndN, preCalPairNMs, overlapMN00);

    //if (!ZOdd) overlapMScheme(pairPMs, basesPM1, ZOdd, P0sP, orbitMProton.size(), overlapMP11);
    //if (!NOdd) overlapMScheme(pairNMs, basesNM1, NOdd, P0sN, orbitMProton.size(), overlapMN11);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "overlap duration: " << duration.count() << " ms" << std::endl;

    start = end;


    Eigen::MatrixXd transformMatrixP0, transformMatrixN0, transformMatrixP1, transformMatrixN1;

    // calculate J-scheme overlap
    transferMatrix(pairPs, pairPMs, A0sP, P0sP, basesJP, basesPM0, transformMatrixP0);
    transferMatrix(pairNs, pairNMs, A0sN, P0sN, basesJN, basesNM0, transformMatrixN0);

    if (!ZOdd) transferMatrix(pairPs, pairPMs, A0sP, P0sP, basesJP, basesPM1, transformMatrixP1);
    if (!NOdd) transferMatrix(pairNs, pairNMs, A0sN, P0sN, basesJN, basesNM1, transformMatrixN1);

    removeUselessBasesJ(basesJP, transformMatrixP0, transformMatrixP1);
    removeUselessBasesJ(basesJN, transformMatrixN0, transformMatrixN1);

    buildBlockJ(basesJP, blockJP);
    buildBlockJ(basesJN, blockJN);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "transfer matrix duration: " << duration.count() << " ms" << std::endl;

    start = end;
    
    matrixElementMtoJ(overlapMP00, transformMatrixP0, transformMatrixP0, overlapJP);
    matrixElementMtoJ(overlapMN00, transformMatrixN0, transformMatrixN0, overlapJN);

    gramSchmidt(overlapJP, blockJP, basesJP, orthogonalBasisP, transformMatrixP0, transformMatrixP1);
    gramSchmidt(overlapJN, blockJN, basesJN, orthogonalBasisN, transformMatrixN0, transformMatrixN1);

    buildBlockJ(basesJP, blockJP);
    buildBlockJ(basesJN, blockJN);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "M to J and gram duration: " << duration.count() << " ms" << std::endl;

    start = end;


    // calculate density matrix element
    Eigen::SparseMatrix<double> oaby6MatrixP, oaby6MatrixN;
    Eigen::MatrixXd qabMatrixP, qabMatrixN;

    oaby6(orbitMProton.size(), orbitP, orbitMProton,
        orbitP, orbitMProton, tbmePP, oaby6MatrixP);
    oaby6(orbitMNeutron.size(), orbitN, orbitMNeutron,
        orbitN, orbitMNeutron, tbmeNN, oaby6MatrixN);

    qab(orbitMProton.size(), orbitMProton, orbitP, qabMatrixP);
    qab(orbitMNeutron.size(), orbitMNeutron, orbitN, qabMatrixN);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "oq duration: " << duration.count() << " ms" << std::endl;

    start = end;


    vector<vector<Eigen::MatrixXd>> fabPM00, fabNM00, fabPM10, fabNM10;
    Eigen::MatrixXd hamMatrixMPP, hamMatrixMNN;
    Eigen::MatrixXd hamMatrixJPP, hamMatrixJNN;
    vector<pair<int, Eigen::MatrixXd>> SRJBlocksP, SRJBlocksN;

    initializeTwoDimVectors(basesPM0.size(), basesPM0.size(), fabPM00);
    initializeTwoDimVectors(basesNM0.size(), basesNM0.size(), fabNM00);
    initializeTwoDimVectors(basesPM0.size(), basesPM1.size(), fabPM10);
    initializeTwoDimVectors(basesNM0.size(), basesNM1.size(), fabNM10);

    generateFabMap(pairPMs, P0sP, basesPM0, basesPM0, orbitMProton.size(), levelEndP, preCalPairPMs, fabPM00);
    generateFabMap(pairNMs, P0sN, basesNM0, basesNM0, orbitMNeutron.size(), levelEndN, preCalPairNMs, fabNM00);

    if (!ZOdd) generateFabMap(pairPMs, P0sP, basesPM0, basesPM1, orbitMProton.size(), levelEndP, preCalPairPMs, fabPM10);
    if (!NOdd) generateFabMap(pairNMs, P0sN, basesNM0, basesNM1, orbitMNeutron.size(), levelEndN, preCalPairNMs, fabNM10);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "single body matrix element duration: " << duration.count() << " ms" << std::endl;

    start = end;

    calHamiltonianMatrix(pairPMs, P0sP, basesPM0, orbitMProton.size(), levelEndP, preCalPairPMs, fabPM00,
        oaby6MatrixP, qabMatrixP, hamMatrixMPP);
    calHamiltonianMatrix(pairNMs, P0sN, basesNM0, orbitMNeutron.size(), levelEndN, preCalPairNMs, fabNM00,
        oaby6MatrixN, qabMatrixN, hamMatrixMNN);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "tbme duration: " << duration.count() << " ms" << std::endl;

    start = end;

    matrixElementMtoJ(hamMatrixMPP, transformMatrixP0, transformMatrixP0, hamMatrixJPP);
    matrixElementMtoJ(hamMatrixMNN, transformMatrixN0, transformMatrixN0, hamMatrixJNN);

    hamMatrixJPP = orthogonalBasisP * hamMatrixJPP * orthogonalBasisP.transpose();
    hamMatrixJNN = orthogonalBasisN * hamMatrixJNN * orthogonalBasisN.transpose();

    Eigen::MatrixXd reducedJP = hamMatrixJPP;
    Eigen::MatrixXd reducedJN = hamMatrixJNN;

    separateReducedJ(reducedJP, basesJP, blockJP, SRJBlocksP);
    separateReducedJ(reducedJN, basesJN, blockJN, SRJBlocksN);

    /*vector<double> eigenValues;
    vector<vector<double>> eigenVectors;

    lanczos(SRJBlocksP[0].second, eigenValues, eigenVectors);

    cout << SRJBlocksP[0].first << "   " <<  eigenValues[0] << "   " <<  eigenValues[1] << endl;

    lanczos(SRJBlocksP[1].second, eigenValues, eigenVectors);

    cout << SRJBlocksP[1].first << "   " <<  eigenValues[0] << "   " <<  eigenValues[1] << endl;*/

    vector<int> ts = {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24};

    vector<vector<Eigen::MatrixXd>> qPi, qNu;
    vector<Eigen::MatrixXd> VIt;
    vector<int> dims;
    dims.push_back(orbitPNumber);
    dims.push_back(orbitPNumber);
    dims.push_back(orbitNNumber);
    dims.push_back(orbitNNumber);
    vector<double> strengthVec;
    strengthVec.push_back(1.0);

    getSVDResult(tbmePN, orbitP, orbitN, dims, ts, qPi, VIt, qNu, strengthVec);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "SVD duration: " << duration.count() << " ms" << std::endl;

    start = end;

    auto rank = qPi[0].size();

    vector<vector<Eigen::MatrixXd>> qabtMMatrixP, qabtMMatrixN;
    vector<vector<Eigen::MatrixXd>> qMMatrixP00, qMMatrixN00, qJMatrixP00, qJMatrixN00;
    vector<vector<Eigen::MatrixXd>> qMMatrixP01, qMMatrixN01, qJMatrixP01, qJMatrixN01;
    vector<vector<Eigen::MatrixXd>> reducedJMatrixP, reducedJMatrixN;

    initializeTwoDimVectors(ts.size(), rank, qabtMMatrixP);
    initializeTwoDimVectors(ts.size(), rank, qabtMMatrixN);
    initializeTwoDimVectors(ts.size(), rank, qMMatrixP00);
    initializeTwoDimVectors(ts.size(), rank, qMMatrixN00);
    initializeTwoDimVectors(ts.size(), rank, qJMatrixP00);
    initializeTwoDimVectors(ts.size(), rank, qJMatrixN00);
    initializeTwoDimVectors(ts.size(), rank, qMMatrixP01);
    initializeTwoDimVectors(ts.size(), rank, qMMatrixN01);
    initializeTwoDimVectors(ts.size(), rank, qJMatrixP01);
    initializeTwoDimVectors(ts.size(), rank, qJMatrixN01);
    initializeTwoDimVectors(ts.size(), rank, reducedJMatrixP);
    initializeTwoDimVectors(ts.size(), rank, reducedJMatrixN);

    qabtJtoM(orbitMProton, orbitP, ts, rank, qPi, qabtMMatrixP);
    qabtJtoM(orbitMNeutron, orbitN, ts, rank, qNu, qabtMMatrixN);

    generateQMMatrix(basesPM0, basesPM0, qabtMMatrixP, fabPM00, qMMatrixP00);
    generateQMMatrix(basesNM0, basesNM0, qabtMMatrixN, fabNM00, qMMatrixN00);

    if (!ZOdd) generateQMMatrix(basesPM0, basesPM1, qabtMMatrixP, fabPM10, qMMatrixP01);
    if (!NOdd) generateQMMatrix(basesNM0, basesNM1, qabtMMatrixN, fabNM10, qMMatrixN01);

    auto loop1 = qMMatrixP00.size();
    auto loop2 = qMMatrixP00[0].size();

    for (int i = 0; i < loop1; i++) {
        for (int j = 0; j < loop2; j++) {
            matrixElementMtoJ(qMMatrixP00[i][j], transformMatrixP0, transformMatrixP0, qJMatrixP00[i][j]);
            matrixElementMtoJ(qMMatrixN00[i][j], transformMatrixN0, transformMatrixN0, qJMatrixN00[i][j]);

            if (!ZOdd) matrixElementMtoJ(qMMatrixP01[i][j], transformMatrixP0, transformMatrixP1, qJMatrixP01[i][j]);
            if (!NOdd) matrixElementMtoJ(qMMatrixN01[i][j], transformMatrixN0, transformMatrixN1, qJMatrixN01[i][j]);
        }
    }


    generateReducedOrthogonalQJ(ts, rank, basesJP, qJMatrixP00, qJMatrixP01, orthogonalBasisP, reducedJMatrixP);
    generateReducedOrthogonalQJ(ts, rank, basesJN, qJMatrixN00, qJMatrixN01, orthogonalBasisN, reducedJMatrixN);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "prepare Q dot Q duration: " << duration.count() << " ms" << std::endl;

    start = end;

    map<int, pair<vector<double>, vector<vector<double>>>> eigens;
    map<int, vector<vector<int>>> basesJPNMap;

    for (auto J : Js) {

        // 0 Pi, 1 Nu, 2 J, 3 M
        vector<vector<int>> basesJPN;

        generateBasesJPNOne(basesJP, basesJN, J, basesJPN);

        Eigen::MatrixXd hPiNuMatrix;

        generateTwoBodyPN(basesJP, basesJN, basesJPN, basesJPN, ts, VIt, reducedJMatrixN, reducedJMatrixP, hPiNuMatrix);

        generateTwoBodyPN(basesJP, basesJN, basesJPN, basesJPN, reducedJP, reducedJN, hPiNuMatrix);

        vector<double> eigenValues;
        vector<vector<double>> eigenVectors;

        lanczos(hPiNuMatrix, eigenValues, eigenVectors);

        cout << eigenValues[0] << eigenValues[1] << endl;

        eigens[J] = {eigenValues, eigenVectors};
        basesJPNMap[J] = basesJPN;
    }

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "spectrum diagonalization duration: " << duration.count() << " ms" << std::endl;

    // calculate BE2

    /*double alphaValue = 4.33195830385423;

    Eigen::MatrixXd qj1j2JP, qj1j2JN, qj1j2MP, qj1j2MN, QopMP00, QopMN00, QopMP01, QopMN01,
    QopJP00, QopJN00, QopJP01, QopJN01, reducedQopJP, reducedQopJN;

    generateQj1j2J(orbitProton, 4, alphaValue, qj1j2JP);
    generateQj1j2J(orbitNeutron, 4, alphaValue, qj1j2JN);

    qj1j2JtoM(orbitMProton, orbitProton, qj1j2JP, 4, qj1j2MP);
    qj1j2JtoM(orbitMNeutron, orbitNeutron, qj1j2JN, 4, qj1j2MN);

    generateQOperatorM(basesPM0, basesPM0, fabPM00, qj1j2MP, QopMP00);
    generateQOperatorM(basesNM0, basesNM0, fabNM00, qj1j2MN, QopMN00);

    generateQOperatorM(basesPM0, basesPM1, fabPM10, qj1j2MP, QopMP01);
    generateQOperatorM(basesNM0, basesNM1, fabNM10, qj1j2MN, QopMN01);

    matrixElementMtoJ(QopMP00, transformMatrixP0, transformMatrixP0, QopJP00);
    matrixElementMtoJ(QopMN00, transformMatrixN0, transformMatrixN0, QopJN00);

    matrixElementMtoJ(QopMP01, transformMatrixP0, transformMatrixP1, QopJP01);
    matrixElementMtoJ(QopMN01, transformMatrixN0, transformMatrixN1, QopJN01);

    generateReducedOrthogonalQopJ(basesJP, QopJP00, QopJP01, orthogonalBasisP, 4, reducedQopJP);
    generateReducedOrthogonalQopJ(basesJN, QopJN00, QopJN01, orthogonalBasisN, 4, reducedQopJN);

    State initialState = {4, eigens[4].first[0], eigens[4].second[0]};
    State finalState = {0, eigens[0].first[0], eigens[0].second[0]};

    double be2 = calculateBE2(basesJP, basesJN, basesJPNMap[0], basesJPNMap[4], 4, initialState, finalState,
        1.5, -0.5, reducedQopJP, reducedQopJN);

    cout << be2 << endl;*/

    // calculate bm1

    /*Eigen::MatrixXd mj1j2JP, mj1j2JN, mj1j2MP, mj1j2MN, MopMP00, MopMN00, MopMP01, MopMN01,
    MopJP00, MopJN00, MopJP01, MopJN01, reducedMopJP, reducedMopJN;

    generateMj1j2J(orbitProton, 2, alphaValue, 3.91, 1.0, 1.0, mj1j2JP);
    generateMj1j2J(orbitNeutron, 2, alphaValue, -2.68, 0.0, 1.0, mj1j2JN);

    qj1j2JtoM(orbitMProton, orbitProton, mj1j2JP, 2, mj1j2MP);
    qj1j2JtoM(orbitMNeutron, orbitNeutron, mj1j2JN, 2, mj1j2MN);

    generateQOperatorM(basesPM0, basesPM0, fabPM00, mj1j2MP, MopMP00);
    generateQOperatorM(basesNM0, basesNM0, fabNM00, mj1j2MN, MopMN00);

    generateQOperatorM(basesPM0, basesPM1, fabPM10, mj1j2MP, MopMP01);
    generateQOperatorM(basesNM0, basesNM1, fabNM10, mj1j2MN, MopMN01);

    matrixElementMtoJ(MopMP00, transformMatrixP0, transformMatrixP0, MopJP00);
    matrixElementMtoJ(MopMN00, transformMatrixN0, transformMatrixN0, MopJN00);

    matrixElementMtoJ(MopMP01, transformMatrixP0, transformMatrixP1, MopJP01);
    matrixElementMtoJ(MopMN01, transformMatrixN0, transformMatrixN1, MopJN01);

    generateReducedOrthogonalQopJ(basesJP, MopJP00, MopJP01, orthogonalBasisP, 2, reducedMopJP);
    generateReducedOrthogonalQopJ(basesJN, MopJN00, MopJN01, orthogonalBasisN, 2, reducedMopJN);

    State mInitialState = {4, eigens[4].first[0], eigens[4].second[0]};
    State mFinalState = {4, eigens[4].first[0], eigens[4].second[0]};

    double bm1 = calculateBM1(basesJP, basesJN, basesJPNMap[4], basesJPNMap[4], 2, mInitialState, mFinalState, reducedMopJP, reducedMopJN);

    cout << bm1 << endl;*/

    return 0;
}