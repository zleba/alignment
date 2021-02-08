#include "read.h"
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

//Load data from Mille binary file
vector<pair<vector<int>,vector<float>>> readData(string fName)
{
    vector<pair<vector<int>,vector<float>>> data;

    ifstream file(fName, std::ifstream::binary); 
    while(file.good()) {
        int nr, tot, dummy;
        file.read ((char*)&tot,4);
        if(!file.good()) break;
        file.read ((char*)&nr,4);
        nr /= 2;

        vector<float> pars(nr);
        file.read ((char*)pars.data(), 4*nr);

        vector<int> indx(nr);
        file.read ((char*)indx.data(), 4*nr);

        file.read ((char*)&dummy,4);

        assert(dummy == 8*nr+4);

        // index vector and pars vector
        data.push_back(make_pair(indx, pars));
    }
    file.close();
    return data;
}

// from the text file
vector<pair<vector<int>,vector<double>>> readConstraints(string fName)
{
    vector<pair<vector<int>,vector<double>>> data;

    //vector<int> indx;
    //vector<double> vals;

    ifstream file(fName); 
    int i = 0;
    while(1) {
        string line;
        getline(file, line);
        if(!file.good()) break;
        size_t startpos = line.find_first_not_of(" \t");
        if(startpos == std::string::npos) continue;
        string str = line.substr( startpos );

        stringstream Str(str);
        if(str[0] == 'C') {
            string w;
            double val;
            Str >> w >> val;
            cout << "C: " << w <<" "<< val << endl;
            data.push_back({{}, {}});
        }
        else {
            int ind;
            double val;
            Str >> ind >> val;
            cout << "V: " << ind <<" "<< val << endl;
            data.back().first.push_back(ind);
            data.back().second.push_back(val);
        }

        //int ind;
        //float val;
        //file >> a >> val;
        //cout << a <<" "<< val << endl;

    }

    file.close();
    return data;

}


//Data to more human readable format
vector<vector<Hit>> convertData(const vector<pair<vector<int>,vector<float>>> &dataOld)
{
    vector<vector<Hit>> data(dataOld.size());
    for(int k = 0; k < dataOld.size(); ++k) {
        auto &d = dataOld[k];
        int i = 1;
        while(i < d.first.size()) {
            Hit hit;
            assert(d.first[i] == 0);
            hit.res = -d.second[i];
            ++i;
            while(d.first[i] > 0) { //Local parameters
                hit.loc.push_back(d.second[i]);
                hit.locIndx.push_back(d.first[i]);
                assert(d.first[i] < 10);
                ++i;
            }
            assert(d.first[i] == 0);
            hit.sigma = d.second[i];
            ++i;
            
            while(d.first[i] > 0) { //Global parameters
                hit.gl.push_back(d.second[i]);
                hit.glIndx.push_back(d.first[i]);
                assert(d.first[i] >= 10);
                ++i;
                if(i >= d.first.size()) break;
            }
            data[k].push_back(hit);
        }
    }

    return data;
}


arma::mat convertConstraints(vector<pair<vector<int>,vector<double>>> con, const vector<map<int,int>> &resIndx)
{
    auto &glIndx  = resIndx[1];
    arma::mat ConM(glIndx.size(), con.size(), arma::fill::zeros);
    for(int ic = 0; ic < con.size(); ++ic) {
        for(int i = 0; i < con[ic].first.size(); ++i)
            ConM(glIndx.at(con[ic].first[i]), ic) = con[ic].second[i];
    }
    return ConM;
}

