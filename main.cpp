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
    string spsFile, pairFile, scFileP, scFileN;
    int Z, N, interactionNumber, orbitPNumber, orbitNNumber, pairTypes, limitNumber;
    vector<Orbit*> orbitProton, orbitNeutron;
    vector<InteractionFile*> interactionFiles;
    vector<Pair*> pairPs, pairNs;
    vector<PairLimit*> limits;

    // a, b, c, d, J, T, value
    map<TBMEJ, map<pair<int, int>, double>> tbmeJMap, tbmePN, tbmePP, tbmeNN;

    // read input files
    readInput(Z, N, spsFile, pairFile, scFileP, scFileN, interactionNumber, interactionFiles);
    readSps(spsFile, orbitPNumber, orbitNNumber, orbitProton, orbitNeutron);

    readPair(pairFile, pairTypes, pairPs, limitNumber, limits);
    readPair(pairFile, pairTypes, pairNs, limitNumber, limits);

    readStructureCoefficient(scFileP, pairPs, pairTypes, orbitPNumber);
    readStructureCoefficient(scFileN, pairNs, pairTypes, orbitNNumber);

    readInteraction(interactionFiles[0]->filename, orbitProton, orbitNeutron, tbmeJMap);
    readInteraction(interactionFiles[1]->filename, orbitProton, orbitNeutron, tbmePN);

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
    generateBasis(pairPs, Z / 2, basis, basesP, limits);
    generateBasis(pairNs, N / 2, basis, basesN, limits);

    generateJi(pairPs, basesP, JisListP);
    generateJi(pairNs, basesN, JisListN);

    generateBasesJ(JisListP, basesP, basesJP);
    generateBasesJ(JisListN, basesN, basesJN);

    //buildBlockJ(basesJP, blockJP);
    //buildBlockJ(basesJN, blockJN);

    initializeOrbitM(orbitProton, orbitMProton);
    initializeOrbitM(orbitNeutron, orbitMNeutron);

    initializePairM(pairPs, orbitProton, orbitMProton, pairPMs, pairJMMapP);
    initializePairM(pairNs, orbitNeutron, orbitMNeutron, pairNMs, pairJMMapN);

    generateMBases(pairPs, pairPMs, Z / 2, basesP, pairJMMapP, basesPM0, basesPM1);
    generateMBases(pairNs, pairNMs, N / 2, basesN, pairJMMapN, basesNM0, basesNM1);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "generate bases duration: " << duration.count() << " ms" << std::endl;

    start = end;


    Eigen::MatrixXd overlapMP00, overlapJP, orthogonalBasisP, overlapMN00, overlapJN, orthogonalBasisN, overlapMP11,
    overlapMN11;

    // calculate M-scheme overlap
    overlapMScheme(pairPMs, basesPM0, orbitMProton.size(), overlapMP00);
    overlapMScheme(pairNMs, basesNM0, orbitMProton.size(), overlapMN00);

    overlapMScheme(pairPMs, basesPM1, orbitMProton.size(), overlapMP11);
    overlapMScheme(pairNMs, basesNM1, orbitMProton.size(), overlapMN11);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "overlap duration: " << duration.count() << " ms" << std::endl;

    start = end;


    Eigen::MatrixXd transformMatrixP0, transformMatrixN0, transformMatrixP1, transformMatrixN1;

    // calculate J-scheme overlap
    transferMatrix(pairPs, pairPMs, basesJP, basesPM0, transformMatrixP0);
    transferMatrix(pairNs, pairNMs, basesJN, basesNM0, transformMatrixN0);

    transferMatrix(pairPs, pairPMs, basesJP, basesPM1, transformMatrixP1);
    transferMatrix(pairNs, pairNMs, basesJN, basesNM1, transformMatrixN1);

    //removeUselessBasesJ(basesJP, transformMatrixP0, transformMatrixP1);
    //removeUselessBasesJ(basesJN, transformMatrixN0, transformMatrixN1);

    buildBlockJ(basesJP, blockJP);
    buildBlockJ(basesJN, blockJN);
    
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


    // deal with interaction
    interactionAntiSymmetric(tbmeJMap, orbitProton);
    interactionAntiSymmetricPN(tbmePN);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "deal with interaction duration: " << duration.count() << " ms" << std::endl;

    start = end;


    // calculate density matrix element
    Eigen::SparseMatrix<double> oaby6MatrixP, oaby6MatrixN;
    Eigen::MatrixXd qabMatrixP, qabMatrixN;

    oaby6(orbitMProton.size(), orbitProton, orbitMProton,
        orbitProton, orbitMProton, tbmeJMap, oaby6MatrixP);
    oaby6(orbitMNeutron.size(), orbitNeutron, orbitMNeutron,
        orbitNeutron, orbitMNeutron, tbmeJMap, oaby6MatrixN);

    qab(orbitMProton.size(), orbitMProton, orbitProton, qabMatrixP);
    qab(orbitMNeutron.size(), orbitMNeutron, orbitNeutron, qabMatrixN);

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

    generateFabMap(pairPMs, basesPM0, basesPM0, orbitMProton.size(), fabPM00);
    generateFabMap(pairNMs, basesNM0, basesNM0, orbitMNeutron.size(), fabNM00);

    generateFabMap(pairPMs, basesPM0, basesPM1, orbitMProton.size(), fabPM10);
    generateFabMap(pairNMs, basesNM0, basesNM1, orbitMNeutron.size(), fabNM10);

    end = std::chrono::high_resolution_clock::now();

    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "single body matrix element duration: " << duration.count() << " ms" << std::endl;

    start = end;

    calHamiltonianMatrix(pairPMs, basesPM0, orbitMProton.size(), fabPM00,
        oaby6MatrixP, qabMatrixP, hamMatrixMPP);
    calHamiltonianMatrix(pairNMs, basesNM0, orbitMNeutron.size(), fabNM00,
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

    cout << eigenValues[0] << eigenValues[1] << endl;

    lanczos(SRJBlocksP[1].second, eigenValues, eigenVectors);

    cout << eigenValues[0] << eigenValues[1] << endl;*/

    vector<int> Js = {0, 4};
    vector<int> ts = {0, 2, 4, 6, 8, 10};

    vector<vector<Eigen::MatrixXd>> qPi, qNu;
    vector<Eigen::MatrixXd> VIt;
    vector<int> dims;
    dims.push_back(orbitPNumber);
    dims.push_back(orbitPNumber);
    dims.push_back(orbitNNumber);
    dims.push_back(orbitNNumber);
    vector<double> strengthVec;
    strengthVec.push_back(1.0);

    getSVDResult(tbmePN, orbitProton, orbitNeutron, dims, ts, qPi, VIt, qNu, strengthVec);

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

    qabtJtoM(orbitMProton, orbitProton, ts, rank, qPi, qabtMMatrixP);
    qabtJtoM(orbitMNeutron, orbitNeutron, ts, rank, qNu, qabtMMatrixN);

    generateQMMatrix(basesPM0, basesPM0, qabtMMatrixP, fabPM00, qMMatrixP00);
    generateQMMatrix(basesNM0, basesNM0, qabtMMatrixN, fabNM00, qMMatrixN00);

    generateQMMatrix(basesPM0, basesPM1, qabtMMatrixP, fabPM10, qMMatrixP01);
    generateQMMatrix(basesNM0, basesNM1, qabtMMatrixN, fabNM10, qMMatrixN01);

    auto loop1 = qMMatrixP00.size();
    auto loop2 = qMMatrixP00[0].size();

    for (int i = 0; i < loop1; i++) {
        for (int j = 0; j < loop2; j++) {
            matrixElementMtoJ(qMMatrixP00[i][j], transformMatrixP0, transformMatrixP0, qJMatrixP00[i][j]);
            matrixElementMtoJ(qMMatrixN00[i][j], transformMatrixN0, transformMatrixN0, qJMatrixN00[i][j]);

            matrixElementMtoJ(qMMatrixP01[i][j], transformMatrixP0, transformMatrixP1, qJMatrixP01[i][j]);
            matrixElementMtoJ(qMMatrixN01[i][j], transformMatrixN0, transformMatrixN1, qJMatrixN01[i][j]);
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