#include "Header.h"
vector<double> solveWithLU(const vector<vector<double>>& A, const vector<double>& f) {
    int N = A.size();

    // Копируем матрицу A и вектор f
    vector<vector<double>> LU = A;
    vector<double> b = f;

    // Прямой ход (LU-разложение с частичным выбором)
    for (int k = 0; k < N; ++k) {
        // Частичный выбор ведущего элемента
        int max_row = k;
        double max_val = abs(LU[k][k]);
        for (int i = k + 1; i < N; ++i) {
            if (abs(LU[i][k]) > max_val) {
                max_val = abs(LU[i][k]);
                max_row = i;
            }
        }

        // Перестановка строк
        if (max_row != k) {
            swap(LU[k], LU[max_row]);
            swap(b[k], b[max_row]);
        }

        // Исключение
        for (int i = k + 1; i < N; ++i) {
            LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < N; ++j) {
                LU[i][j] -= LU[i][k] * LU[k][j];
            }
        }
    }

    // Решение Ly = b
    vector<double> y(N);
    for (int i = 0; i < N; ++i) {
        y[i] = b[i];
        for (int j = 0; j < i; ++j) {
            y[i] -= LU[i][j] * y[j];
        }
    }

    // Решение Ux = y
    vector<double> x(N);
    for (int i = N - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < N; ++j) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] /= LU[i][i];
    }

    return x;
}