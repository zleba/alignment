#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <armadillo>

using namespace std;

struct Hit {
    float sigma, res;
    vector<int> glIndx, locIndx;
    vector<float> gl, loc;
};


vector<pair<vector<int>,vector<float>>> readData(string fName);
vector<pair<vector<int>,vector<float>>> readConstraints(string fName); //from text file
vector<map<int,int>> getIndexes(const  vector<vector<Hit>>  &data);
vector<vector<Hit>> convertData(const vector<pair<vector<int>,vector<float>>> &dataOld);
vector<vector<Hit>> generateData(int nTracks, int nHits, int nLoc, int nGl);
arma::fmat convertConstraints(vector<pair<vector<int>,vector<float>>> con, const vector<map<int,int>> &resIndx);


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

        data.push_back(make_pair(indx, pars));
    }
    file.close();
    return data;
}

vector<pair<vector<int>,vector<float>>> readConstraints(string fName)
{
    vector<pair<vector<int>,vector<float>>> data;

    //vector<int> indx;
    //vector<float> vals;

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
            float val;
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


//Get range of local and global indexes
vector<map<int,int>> getIndexes(const  vector<vector<Hit>>  &data)
{
    map<int,int> locIndx, glIndx;

    for(auto & dd : data) {
        for(auto &d : dd ) {
            for(auto ind : d.locIndx) locIndx[ind] = 1;
            for(auto ind : d.glIndx)  glIndx[ind] = 1;
        }
    }

    int i = 0;
    for(auto &l : locIndx) l.second = i++;
    int j = 0;
    for(auto &g : glIndx)  g.second = j++;

    return {locIndx, glIndx};
}


float Rand() { return rand()/(RAND_MAX+0.);}

//Generate random data for testing
vector<vector<Hit>> generateData(int nTracks, int nHits, int nLoc, int nGl)
{
    vector<vector<Hit>> data(nTracks);
    for(auto & dd : data) {
        dd.resize(nHits);
        for(auto & d : dd) {
            d.res = Rand();
            d.sigma = 0.10;
            d.gl     = { Rand() };
            d.glIndx = { 10};
            d.loc     = { Rand(), Rand() };
            d.locIndx = { 1, 2};
        }
    }
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
            hit.res = d.second[i];
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

struct Fitter {

    arma::fmat C1, Con;
    arma::fvec g1;

    vector<arma::fmat> Gamma, G;
    vector<arma::fvec> beta;

    vector<vector<Hit>> data;

    map<int,int> locIndx, glIndx;

    //void getChi2();
    float getChi2(arma::fvec rGl =  arma::zeros<arma::fvec>(0) );
    void MoveLocal();
    arma::fvec reduceData(const vector<vector<Hit>> &dataNow, const vector<map<int,int>> &resIndx);

};


arma::fmat convertConstraints(vector<pair<vector<int>,vector<float>>> con, const vector<map<int,int>> &resIndx)
{
    auto &glIndx  = resIndx[1];
    arma::fmat Con(glIndx.size(), con.size(), arma::fill::zeros);
    for(int ic = 0; ic < con.size(); ++ic) {
        for(int i = 0; i < con[ic].first.size(); ++i)
            Con(glIndx.at(con[ic].first[i]), ic) = con[ic].second[i];
    }
    return Con;
}


float Fitter::getChi2(arma::fvec rGl)
{
    int nPt = 0;
    double chi2 = 0;
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        arma::fvec qVec = - Gamma[itr].i() * beta[itr];
        if(rGl.n_rows > 0)
            qVec -= G[itr].t()*rGl;

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            double z = d.res;
            for(int j = 0; j < qVec.n_rows; ++j)
                z += d.loc[j] * qVec[j];
            for(int j = 0; j < d.gl.size(); ++j) {
                int glI = glIndx.at(d.glIndx[j]);
                z += d.gl[j]  * rGl[glI];
            }


            chi2 += pow(z/d.sigma, 2);
            ++nPt;
        }
    }
    //cout << "chi2 = " << chi2 << " "<< chi2/nPt << endl;
    return chi2;
}


void Fitter::MoveLocal()
{
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        arma::fvec qVec = - Gamma[itr].i() * beta[itr];

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            double z = d.res;
            for(int j = 0; j < qVec.n_rows; ++j)
                z += d.loc[j] * qVec[j];

            d.res = z;
        }
    }
}





//Calculate the important matrices
arma::fvec Fitter::reduceData(const vector<vector<Hit>> &dataNow, const vector<map<int,int>> &resIndx)
{
    data = dataNow;

    locIndx = resIndx[0];
    glIndx  = resIndx[1];
    C1.zeros(glIndx.size(), glIndx.size());
    g1.zeros(glIndx.size());

    Gamma.resize(data.size());
    G.resize(data.size());
    beta.resize(data.size());
    
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        Gamma[itr].zeros(locIndx.size(),locIndx.size());
        G[itr].zeros(glIndx.size(),locIndx.size());
        beta[itr].zeros(locIndx.size());

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            float sigmaInv = 1/(d.sigma*d.sigma);

            //Fill C1
            for(int i = 0; i < d.gl.size(); ++i)
            for(int j = 0; j < d.gl.size(); ++j) {
                C1(glIndx.at(d.glIndx[i]), glIndx.at(d.glIndx[j])) += d.gl[i]*d.gl[j] *sigmaInv;
            }

            //Fill Gamma
            for(int i = 0; i < d.loc.size(); ++i)
            for(int j = 0; j < d.loc.size(); ++j) {
                Gamma[itr](locIndx.at(d.locIndx[i]), locIndx.at(d.locIndx[j])) += d.loc[i]*d.loc[j] *sigmaInv;
            }

            //Fill G
            for(int i = 0; i < d.gl.size(); ++i)
            for(int j = 0; j < d.loc.size(); ++j) {
                G[itr](glIndx.at(d.glIndx[i]), locIndx.at(d.locIndx[j])) += d.gl[i]*d.loc[j] *sigmaInv;
            }

            //Fill beta
            for(int i = 0; i < d.loc.size(); ++i) {
                beta[itr](locIndx.at(d.locIndx[i])) += d.loc[i]*d.res *sigmaInv;
            }

            //Fill g1
            for(int i = 0; i < d.gl.size(); ++i) {
                g1(glIndx.at(d.glIndx[i])) += d.gl[i]*d.res *sigmaInv;
            }
        }
    }

    arma::fmat C2(glIndx.size(), glIndx.size(), arma::fill::zeros);
    arma::fvec g2(glIndx.size());
    //Sum the matrices
    for(int itr = 0; itr < data.size(); ++itr) {
        //arma::fmat Start =  - G[itr] * arma::inv_sympd(Gamma[itr]);
        arma::fmat Start =  - G[itr] * Gamma[itr].i();
        //arma::fmat GammaI = Gamma[itr].i();

        C2 += Start * G[itr].t();
        g2 += Start * beta[itr];
    }

    arma::fmat C = C1 + C2;
    arma::fvec g = g1 + g2;

    cout << "cond " << arma::cond(C) << endl;
    //cout << "vals " << C(20,20) <<" "<< g(20) << endl;



    arma::fvec r = - C.i() * g;

    return r;
    //cout << r << endl;
}

void test()
{
    auto data = generateData(6, 20, 2, 3);
    auto indx = getIndexes(data);

    Fitter fitter;
    auto r = fitter.reduceData(data, indx);
    //fitter.MoveLocal();
    //r = fitter.reduceData(data, indx);

    for(float a = -0.001; a < 0.001; a +=0.00001) {
        arma::fvec v({a});
        float chi = fitter.getChi2(v);
        cout << "a: " << a << " "<< chi << endl;
        
    }

}

int main()
{
    test();
    return 0;

    auto con = readConstraints("/home/radek/Downloads/pede/target/mp2con.txt"); //from text file
    return 0;

    auto dataOrg = readData("/home/radek/Downloads/pede/target/mp2tst.bin");



    auto data =  convertData(dataOrg);
    auto indx = getIndexes(data);
    //auto indxO = getIndexes(dataOrg);
    //assert(indxO == indx);

    Fitter fitter;
    fitter.Con = convertConstraints(con, indx);
    fitter.reduceData(data, indx);
    fitter.getChi2();

    fitter.MoveLocal();
    auto r = fitter.reduceData(data, indx);
    cout << "r20 = " << r(20) << endl;
    fitter.getChi2(r);

    return 0;
}
