#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <vector>
#include <map>
#include <armadillo>

#include "read.h"

using namespace std;



vector<map<int,int>> getIndexes(const  vector<vector<Hit>>  &data);
vector<vector<Hit>> generateData(int nTracks, int nHits, int nLoc, int nGl);
arma::mat convertConstraints(vector<pair<vector<int>,vector<double>>> con, const vector<map<int,int>> &resIndx);




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


double Rand() { return rand()/(RAND_MAX+0.);}


double Gauss()
{
    double s = 0;
    for(int i = 0; i < 12; ++i)
        s += Rand()-0.5;
    return s;
}

//Generate random data for testing
vector<vector<Hit>> generateDataToy(int nTracks, vector<double> parsIn)
{
    vector<vector<Hit>> data(nTracks);
    for(auto & dd : data) {
        dd.resize(parsIn.size());

        double s = 2*Rand() -1; //-1, 1
        double a = 5*(Rand() - 0.5);

        vector<double> v(parsIn.size());
        for(int i = 0; i < v.size(); ++i)
            v[i] = a + i*s + Gauss() + parsIn[i];

        //very simple tracking
        double sN = (v.back() - v.front()) / (v.size()-1);
        double aN = v.front();


        for(int i = 0; i < v.size(); ++i) {
            dd[i].res = (aN+i*sN) - v[i];
            dd[i].sigma = 1;
            dd[i].gl     = {1.};
            dd[i].glIndx = {10+i};
            dd[i].loc     = { 1., 1.*i };
            dd[i].locIndx = { 1, 2};
        }
    }
    return data;
}





struct Fitter {

    arma::mat C1, ConM;
    arma::vec g1, ConR;

    vector<arma::mat> Gamma, G;
    vector<arma::vec> beta;

    vector<vector<Hit>> data;

    map<int,int> locIndx, glIndx;

    double getChi2(arma::vec rGl);

    void MoveLocalGlobal(arma::vec rGl=arma::vec());
    void MoveGlobal(arma::vec rGl);
    arma::vec reduceData(const vector<vector<Hit>> &dataNow, const vector<map<int,int>> &resIndx);

};



double Fitter::getChi2(arma::vec rGl)
{
    int nPt = 0;
    double chi2 = 0;
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        arma::vec qVec = - Gamma[itr].i() * beta[itr];

        //cout << "qVecNorm : " << arma::norm(qVec) << endl;

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            //loop over local parameters
            double z = +d.res;
            for(int j = 0; j < d.loc.size(); ++j) {
                int lI = locIndx.at(d.locIndx[j]);
                z += d.loc[j] * qVec[lI];
            }

            //loop over global parameters
            for(int j = 0; j < d.gl.size(); ++j) {
                int glI = glIndx.at(d.glIndx[j]);
                z += d.gl[j]  * rGl[glI];
            }



            chi2 += pow(z/d.sigma, 2);
            ++nPt;
        }

        nPt -= dd[0].loc.size(); //subtract local degrees of freedom
    }
    cout << "chi2 = " << chi2 << " "<< chi2/nPt << endl;
    return chi2;
}


void Fitter::MoveLocalGlobal(arma::vec rGl)
{
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        arma::vec qVec = - Gamma[itr].i() * beta[itr];

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            double z = d.res;
            for(int j = 0; j < d.loc.size(); ++j) {
                int lI = locIndx.at(d.locIndx[j]);
                z += d.loc[j] * qVec[lI];
            }

            if(rGl.n_rows > 0) {
                //loop over global parameters
                for(int j = 0; j < d.gl.size(); ++j) {
                    int glI = glIndx.at(d.glIndx[j]);
                    z += d.gl[j]  * rGl[glI];
                }
            }


            d.res = z;
        }
    }
}


void Fitter::MoveGlobal(arma::vec rGl)
{
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            double z = d.res;
            //loop over global parameters
            for(int j = 0; j < d.gl.size(); ++j) {
                int glI = glIndx.at(d.glIndx[j]);
                z += d.gl[j]  * rGl[glI];
            }

            d.res = z;
        }
    }
}







//Calculate the important matrices
arma::vec Fitter::reduceData(const vector<vector<Hit>> &dataNow, const vector<map<int,int>> &resIndx)
{
    data = dataNow;

    locIndx = resIndx[0];
    glIndx  = resIndx[1];
    C1.zeros(glIndx.size(), glIndx.size());
    g1.zeros(glIndx.size());

    Gamma.resize(data.size());
    G.resize(data.size());
    beta.resize(data.size());
    
    // event loop
    for(int itr = 0; itr < data.size(); ++itr) {
        auto & dd = data[itr];

        Gamma[itr].zeros(locIndx.size(),locIndx.size());
        G[itr].zeros(glIndx.size(),locIndx.size());
        beta[itr].zeros(locIndx.size());

        for(int ihit = 0; ihit < dd.size(); ++ihit) {
            auto &d = dd[ihit];

            double sigmaInv = 1/(d.sigma*d.sigma);

            //Fill C1 (stable)
            for(int i = 0; i < d.gl.size(); ++i)
            for(int j = 0; j < d.gl.size(); ++j) {
                C1(glIndx.at(d.glIndx[i]), glIndx.at(d.glIndx[j])) += d.gl[i]*d.gl[j] *sigmaInv;
            }

            //Fill Gamma (stable)
            for(int i = 0; i < d.loc.size(); ++i)
            for(int j = 0; j < d.loc.size(); ++j) {
                Gamma[itr](locIndx.at(d.locIndx[i]), locIndx.at(d.locIndx[j])) += d.loc[i]*d.loc[j] *sigmaInv;
            }

            //Fill G (stable)
            for(int i = 0; i < d.gl.size(); ++i)
            for(int j = 0; j < d.loc.size(); ++j) {
                G[itr](glIndx.at(d.glIndx[i]), locIndx.at(d.locIndx[j])) += d.gl[i]*d.loc[j] *sigmaInv;
            }

            //Fill beta (res-dep)
            for(int i = 0; i < d.loc.size(); ++i) {
                beta[itr](locIndx.at(d.locIndx[i])) += d.loc[i]*d.res *sigmaInv;
            }

            //Fill g1 (res-dep)
            for(int i = 0; i < d.gl.size(); ++i) {
                g1(glIndx.at(d.glIndx[i])) += d.gl[i]*d.res *sigmaInv;
            }
        }
    }

    arma::mat C2(glIndx.size(), glIndx.size(), arma::fill::zeros);
    arma::vec g2(glIndx.size());

    //Sum the matrices
    for(int itr = 0; itr < data.size(); ++itr) {
        //arma::mat Start =  - G[itr] * arma::inv_sympd(Gamma[itr]);
        arma::mat Start =  - G[itr] * Gamma[itr].i();
        //arma::mat GammaI = Gamma[itr].i();

        C2 += Start * G[itr].t();
        g2 += Start * beta[itr];
    }

    arma::mat C = C1 + C2;
    arma::vec g = g1 + g2;

    cout << "cond " << arma::cond(C) << endl;
    cout << "Constraint " << ConM.n_cols <<" "<< ConM.n_rows << endl;

    arma::mat Ce(C.n_rows + ConM.n_cols, C.n_rows + ConM.n_cols,  arma::fill::zeros);
    cout << "Size1 " << ConM.n_cols << endl;
    Ce(arma::span(0,C.n_rows-1), arma::span(0,C.n_cols-1)) = C;

    if(ConM.n_cols > 0) {
        cout << "Size2 " << ConM.n_cols << endl;
        Ce(arma::span(0,C.n_rows-1), arma::span(C.n_cols,Ce.n_cols-1)) = ConM;
        cout << "Size3 " << ConM.n_cols << endl;
        Ce(arma::span(C.n_rows,Ce.n_rows-1), arma::span(0,C.n_cols-1)) = ConM.t();
        cout << "Size4 " << ConM.n_cols << endl;
    }
    cout << "Rank " << C.n_rows << " "<< Ce.n_rows <<" : " <<  arma::rank(C) <<" "<<  arma::rank(Ce) << endl;

    arma::vec ge(g.n_rows + ConM.n_cols, arma::fill::zeros);
    cout << "Radecek1 " << endl;
    ge(arma::span(0,g.n_rows-1)) = g;
    ge(arma::span(g.n_rows, ge.n_rows-1)) = -ConR;
    cout << "Radecek2 " << endl;


    //cout << ConM << endl; 
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, Ce);
    cout << "Helenka " << eigval << endl;


    //arma::vec r = - C.i()  * g;
    arma::vec r = - Ce.i() * ge;
    cout << "Radecek3 " << endl;
    cout << "Done" << endl;

    return r(arma::span(0, g.n_rows-1));
    //cout << r << endl;
}

void test()
{
    //auto data = generateDataToy(100000, {-0.4, 0.3, 0.6,-0.5});
    auto data = generateDataToy(100000, {0.3, -0.4, 0.5});
    auto indx = getIndexes(data);

    cout << "Reading done " << indx[0].size() <<" "<< indx[1].size() << endl;
    Fitter fitter;




    fitter.ConM = convertConstraints( {
            //make_pair(vector<int>({10,11,12,13}),vector<double>({1.,1.,1.,1.})),
            //make_pair(vector<int>({10,11,12,13}),vector<double>({0.,1.,2.,3.}))

            make_pair(vector<int>({10,11,12}),vector<double>({1.,1.,1.})),
            make_pair(vector<int>({10,11,12}),vector<double>({0.,1.,2.}))

            }, indx);
    fitter.ConR = arma::vec({+0.4,+0.6});

    auto r = fitter.reduceData(data, indx);
    cout << r << endl;


    double chi2 = fitter.getChi2(r);

    //second deriv
    arma::vec rd=0*r;
    rd[0] = 0.00001;
    rd[1] =-0.00002;
    rd[2] = 0.00001;
    double diff1 = (fitter.getChi2(r+rd) - fitter.getChi2(r)) / 0.00001;
    double diff2 = (fitter.getChi2(r-rd) - fitter.getChi2(r)) / 0.00001;
    cout << "diff " << diff1 <<" "<< diff2 << endl;


    cout << "hela " << chi2 << endl;

    //fitter.MoveLocalGlobal(r);
    fitter.MoveGlobal(r);
    fitter.ConR = arma::vec({+0.0,+0.0});
    auto rr = fitter.reduceData(fitter.data, indx);
    cout << rr << endl;
    chi2 = fitter.getChi2(rr);

}

void testFile()
{
    auto data = convertData(readData("/home/radek/Downloads/pede/pede/mp2tst.bin"));
    auto con  = readConstraints("/home/radek/Downloads/pede/pede/mp2con.txt"); //from text file
    auto indx = getIndexes(data);

    cout << "Reading done " << indx[0].size() <<" "<< indx[1].size() << endl;
    Fitter fitter;

    fitter.ConM = convertConstraints( con, indx);
    fitter.ConR = arma::vec({+0.0,+0.0});

    auto r = fitter.reduceData(data, indx);
    cout << r << endl;


    double chi2 = fitter.getChi2(r);

    /*
    //second deriv
    arma::vec rd=0*r;
    rd[0] = 0.00001;
    rd[1] =-0.00002;
    rd[2] = 0.00001;
    double diff1 = (fitter.getChi2(r+rd) - fitter.getChi2(r)) / 0.00001;
    double diff2 = (fitter.getChi2(r-rd) - fitter.getChi2(r)) / 0.00001;
    cout << "diff " << diff1 <<" "<< diff2 << endl;
    */


    cout << "hela " << chi2 << endl;

    fitter.MoveLocalGlobal(r);
    //fitter.MoveGlobal(r);
    fitter.ConR = arma::vec({+0.0,+0.0});
    auto rr = fitter.reduceData(fitter.data, indx);
    cout << rr << endl;
    chi2 = fitter.getChi2(rr);

}






int main()
{
    testFile();
    return 0;
    test();
    return 0;

    auto con = readConstraints("/home/radek/Downloads/pede/pede/mp2con.txt"); //from text file
    //return 0;

    auto dataOrg = readData("/home/radek/Downloads/pede/pede/mp2tst.bin");



    auto data =  convertData(dataOrg);
    auto indx = getIndexes(data);
    //auto indxO = getIndexes(dataOrg);
    //assert(indxO == indx);

    Fitter fitter;
    fitter.ConM = convertConstraints(con, indx);
    arma::vec r = fitter.reduceData(data, indx);
    double chi2 = fitter.getChi2(r);
    cout << "chi2 " << chi2 << endl;
    //return 0;

    fitter.MoveLocalGlobal();
    auto rr = fitter.reduceData(data, indx);
    cout << "r20 = " << rr(20) << endl;
    cout << "chi2=" << fitter.getChi2(rr) << endl;

    return 0;
}
