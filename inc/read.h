#pragma once

#include <vector>
#include <map>
#include <utility>
#include <armadillo>


struct Hit {
    double sigma, res;
    std::vector<int> glIndx, locIndx;
    std::vector<double> gl, loc;
};

std::vector<std::pair<std::vector<int>,std::vector<float>>> readData(std::string fName); //from Mille binary file

std::vector<std::pair<std::vector<int>,std::vector<double>>> readConstraints(std::string fName); //from text file

std::vector<std::vector<Hit>> convertData(const std::vector<std::pair<std::vector<int>,std::vector<float>>> &dataOld);

arma::mat convertConstraints(std::vector<std::pair<std::vector<int>,std::vector<double>>> con, const std::vector<std::map<int,int>> &resIndx);
