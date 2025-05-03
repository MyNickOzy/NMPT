#include "Header.h"
vector<double> createVectorF(const vector<vector<double>>& A) {
    int N = A.size();
    vector<double> f(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            f[i] += A[i][j]; // Умножение на 1 (x*[j] = 1)
        }
    }
    return f;
}