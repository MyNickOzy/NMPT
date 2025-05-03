#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <sstream>

using namespace std;
using namespace std::chrono;

vector<vector<double>> createMatrixA(int N);
vector<double> createVectorF(const vector<vector<double>>& A);
double euclideanNorm(const vector<double>& v);
double calculateError(const vector<double>& x);
vector<double> solveWithLU(const vector<vector<double>>& A, const vector<double>& f);
vector<double> solveWithQR(const vector<vector<double>>& A, const vector<double>& f);
