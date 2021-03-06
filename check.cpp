#include <iostream>
#include <armadillo>
#include <chrono>

using namespace std;

extern "C" void sqminv_(double *v, double *b, int *n, int *nrank, double *diag, int *next);
extern "C" void sqminl_(double *v, double *b, int *n, int *nrank, double *diag, int *next);

vector<double> toVector(arma::mat test)
{
    int n = test.n_rows;
    vector<double> testV((n*(n+1))/2);
    int k = 0;
    for(int j = 0; j < n; ++j)
        for(int i = 0; i <= j; ++i)
            testV[k++] = test(i, j);
    return testV;
}

arma::mat toMatrix(vector<double> testV)
{
    int n = sqrt(testV.size()*2);
    //vector<double> testV((n*(n+1))/2);
    arma::mat test(n,n);
    int k = 0;
    for(int j = 0; j < n; ++j)
        for(int i = 0; i <= j; ++i) {
            test(i,j) = test(j,i) = testV[k++];
            }
    return test;
}


int getTime(int n)
{
    //const int n = 1024*2*2;
    arma::mat test;
    test.randn(n, n);
    test *= test.t();


    vector<double> testV = toVector(test);

    //double b[n], diag[n];
    double *diag = new double[n];
    double *b    = new double[n];
    int *next = new int[n];

    int nrank, nNow = n;

    //cout << "Init done " << testV[5]<< endl;
    auto begin = std::chrono::steady_clock::now();
    arma::mat C =  arma::inv_sympd(test);
    //arma::mat C =  test.i();

    //sqminl_(testV.data(), b, &nNow, &nrank, diag, next);

    auto end = std::chrono::steady_clock::now();

    delete [] diag;
    delete [] b;
    delete [] next;

    return std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
}





int main(int argc, char **argv)
{
    for(int n = 128; n <= 18000; n *= 2) {
        int t =  getTime(n);
        cout << n <<" : " << t << endl;
    //std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms]" << std::endl;
    }

    //arma::mat C =  test.i();
    //return 0;
    //return 0;
    
    //arma::mat r = toMatrix(testV);
    //cout << "Printing" << endl;
    //cout << (r*test)  << endl;

    return 0;
}
