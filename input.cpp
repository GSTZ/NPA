#include "input.h"


using namespace std;

int readSps(const string& filename, int& orbitProtonNumber, int& orbitNeutronNumber,
    vector<Orbit*>& orbitProton, vector<Orbit*>& orbitNeutron) {
    ifstream infile(filename);
    if (infile.is_open()) {
        int i;
        string line, tem;
        getline(infile, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        istringstream iss(line);
        iss >> tem;
        orbitProtonNumber = stoi(tem);
        iss >> tem;
        orbitNeutronNumber = stoi(tem);
        iss.clear();

        for (i = 0; i < orbitProtonNumber; i++) {
            int n, l, j;
            getline(infile, line);
            if (!line.empty() && line.back() == '\r') line.pop_back();
            iss.str(line);
            iss >> n >> l >> j;
            auto orbit = new Orbit(n, l, j);
            orbitProton.push_back(orbit);
            iss.clear();
        }

        for (i = 0; i < orbitNeutronNumber; i++) {
            int n, l, j;
            getline(infile, line);
            if (!line.empty() && line.back() == '\r') line.pop_back();
            iss.str(line);
            iss >> n >> l >> j;
            auto orbit = new Orbit(n, l, j);
            orbitNeutron.push_back(orbit);
            iss.clear();
        }

        infile.close();
    }

    return 0;
}


int readInput(int& Z, int& N, string& spsFile, string& pairFile, int& interactionNumber,
    vector<InteractionFile*>& interactionFiles) {

    ifstream infile("shell.dat");
    if (infile.is_open()) {
        int i;
        string line, tem;

        getline(infile, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        istringstream iss(line);
        iss >> tem;
        Z = stoi(tem);
        iss >> tem;
        N = stoi(tem);
        iss.clear();

        getline(infile, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        iss.str(line);
        iss >> spsFile;
        iss.clear();

        getline(infile, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        iss.str(line);
        iss >> pairFile;
        iss.clear();

        getline(infile, line);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        iss.str(line);
        iss >> tem;
        interactionNumber = stoi(tem);
        iss.clear();


        for (i = 0; i < interactionNumber; i++) {
            string filename, core, A, scale;
            getline(infile, line);
            if (!line.empty() && line.back() == '\r') line.pop_back();
            iss.str(line);
            iss >> filename >> core >> A >> scale;
            auto infi = new InteractionFile(filename, stoi(core), stoi(A), stod(scale));
            interactionFiles.push_back(infi);
            iss.clear();
        }
    }

    return 0;
}


int readPair(const string& filename, int& pairTypes, vector<Pair*>& pairs, int& limitNumber,
    vector<PairLimit*>& limits) {
    ifstream infile(filename);
    if (infile.is_open()) {
        int i;
        string line1, line2, tem;
        getline(infile, line1);
        if (!line1.empty() && line1.back() == '\r') line1.pop_back();
        istringstream iss1(line1), iss2;
        iss1 >> tem;
        pairTypes = stoi(tem);
        iss1.clear();

        getline(infile, line1);
        if (!line1.empty() && line1.back() == '\r') line1.pop_back();
        getline(infile, line2);
        if (!line2.empty() && line2.back() == '\r') line2.pop_back();
        iss1.str(line1);
        iss2.str(line2);
        for (i = 0; i < pairTypes; i++) {
            string j, parity;
            iss1 >> j;
            iss2 >> parity;
            auto pair = new Pair();
            pair->index = i;
            pair->j = stoi(j);
            if (parity == "1") pair->parity = Parity::Positive;
            else pair->parity = Parity::Negative;
            pairs.push_back(pair);
        }
        iss1.clear();
        iss2.clear();

        getline(infile, line1);
        if (!line1.empty() && line1.back() == '\r') line1.pop_back();
        iss1.str(line1);
        iss1 >> tem;
        limitNumber = stoi(tem);
        iss1.clear();

        for (i = 0; i < limitNumber; i++) {
            int limitMin, limitMax;
            getline(infile, line1);
            if (!line1.empty() && line1.back() == '\r') line1.pop_back();
            iss1.str(line1);
            iss1 >> tem;
            limitMin = stoi(tem);
            iss1 >> tem;
            limitMax = stoi(tem);
            vector<int> indexes;
            while (iss1 >> tem) {
                indexes.push_back(stoi(tem)  - 1); // in C++, index start from 0
            }
            iss1.clear();
            auto limit = new PairLimit(limitMin, limitMax, indexes);
            limits.push_back(limit);
        }
    }
    return 0;
}