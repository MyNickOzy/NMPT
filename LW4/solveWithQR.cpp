#include "Header.h"
vector<double> solveWithQR(const vector<vector<double>>& A, const vector<double>& f) {
    int N = A.size();

    // Инициализация матриц Q и R
    vector<vector<double>> Q(N, vector<double>(N, 0.0));
    vector<vector<double>> R(N, vector<double>(N, 0.0));

    // QR-разложение методом Грама-Шмидта
    for (int j = 0; j < N; ++j) {
        // Копируем j-й столбец A в v
        vector<double> v(N);
        for (int i = 0; i < N; ++i) {
            v[i] = A[i][j];
        }

        // Ортогонализация относительно предыдущих столбцов Q
        for (int k = 0; k < j; ++k) {
            // Вычисляем R[k][j] как скалярное произведение Q[:,k] и A[:,j]
            R[k][j] = 0.0;
            for (int i = 0; i < N; ++i) {
                R[k][j] += Q[i][k] * A[i][j];
            }

            // Вычитаем проекцию
            for (int i = 0; i < N; ++i) {
                v[i] -= R[k][j] * Q[i][k];
            }
        }

        // Нормализация
        R[j][j] = euclideanNorm(v);
        for (int i = 0; i < N; ++i) {
            Q[i][j] = v[i] / R[j][j];
        }
    }

    // Решение системы Rx = Q^T b
    vector<double> QTb(N, 0.0);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            QTb[i] += Q[j][i] * f[j];
        }
    }

    // Обратная подстановка
    vector<double> x(N);
    for (int i = N - 1; i >= 0; --i) {
        x[i] = QTb[i];
        for (int j = i + 1; j < N; ++j) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}