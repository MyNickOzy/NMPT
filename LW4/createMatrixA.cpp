#include "Header.h"
vector<vector<double>> createMatrixA(int N) {
    vector<vector<double>> A(N, vector<double>(N));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                A[i][j] = 100.0;
            }
            else {
                A[i][j] = 10.0 + 0.1 * (i + 1) - 0.2 * (j + 1);
            }
        }
    }
    return A;
}