#include <iostream>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <chrono>
#include <numeric>

using namespace std::chrono;
using namespace std;

// Вычисление нормы вектора
double vectorNorm(const vector<double>&v) {
    double sum = 0.0;
    for (double val : v) {
        sum += val * val;
    }
    return sqrt(sum);
}

// Вычисление погрешности p = ||x - x*|| / ||x*||
double calculateError(const vector<double>& x, const vector<double>& x_exact) {
    vector<double> diff(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        diff[i] = x[i] - x_exact[i];
    }
    return vectorNorm(diff) / vectorNorm(x_exact);
}

// LU-разложение матрицы A (Doolittle algorithm)
void luDecomposition(const vector<vector<double>>& A, vector<vector<double>>& L, vector<vector<double>>& U) {
    int n = A.size();
    L = vector<vector<double>>(n, vector<double>(n, 0.0));
    U = vector<vector<double>>(n, vector<double>(n, 0.0));

    for (int i = 0; i < n; ++i) {
        // Верхняя треугольная матрица U
        for (int k = i; k < n; ++k) {
            double sum = 0.0;
            for (int j = 0; j < i; ++j) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - sum;
        }

        // Нижняя треугольная матрица L
        for (int k = i; k < n; ++k) {
            if (i == k) {
                L[i][i] = 1.0; // Диагональ L заполняется единицами
            }
            else {
                double sum = 0.0;
                for (int j = 0; j < i; ++j) {
                    sum += L[k][j] * U[j][i];
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Решение Ly = f (прямая подстановка)
vector<double> forwardSubstitution(const vector<vector<double>>& L, const vector<double>& f) {
    int n = L.size();
    vector<double> y(n, 0.0);
    for (int i = 0; i < n; ++i) {
        y[i] = f[i];
        for (int j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i]; // L[i][i] = 1, можно опустить
    }
    return y;
}

// Решение Ux = y (обратная подстановка)
vector<double> backwardSubstitution(const vector<vector<double>>& U, const vector<double>& y) {
    int n = U.size();
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
    return x;
}

// Решение Ax = f методом LU-разложения
vector<double> solveLU(const vector<vector<double>>& A, const vector<double>& f) {
    int n = A.size();
    vector<vector<double>> L, U;
    luDecomposition(A, L, U);
    vector<double> y = forwardSubstitution(L, f);
    vector<double> x = backwardSubstitution(U, y);
    return x;
}
// QR-разложение методом Грама-Шмидта
void qrDecomposition(const vector<vector<double>>& A, vector<vector<double>>& Q, vector<vector<double>>& R) {
    int n = A.size();
    Q = vector<vector<double>>(n, vector<double>(n, 0.0));
    R = vector<vector<double>>(n, vector<double>(n, 0.0));

    vector<vector<double>> U = A; // Копия матрицы A

    for (int j = 0; j < n; ++j) {
        // Ортогонализация столбцов
        vector<double> v = U[j];
        for (int i = 0; i < j; ++i) {
            // Проекция v на Q[i]
            R[i][j] = 0.0;
            for (int k = 0; k < n; ++k) {
                R[i][j] += Q[k][i] * U[j][k];
            }
            // Вычитание проекции
            for (int k = 0; k < n; ++k) {
                v[k] -= R[i][j] * Q[k][i];
            }
        }

        // Нормализация
        double norm = 0.0;
        for (double val : v) {
            norm += val * val;
        }
        norm = sqrt(norm);
        R[j][j] = norm;

        // Заполнение Q
        for (int k = 0; k < n; ++k) {
            Q[k][j] = v[k] / norm;
        }
    }
}

// Решение системы Rx = Q^T f (обратная подстановка)
vector<double> solveQR(const vector<vector<double>>& A, const vector<double>& f) {
    int n = A.size();
    vector<vector<double>> Q, R;
    qrDecomposition(A, Q, R);

    // Вычисляем Q^T * f
    vector<double> qt_f(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            qt_f[i] += Q[j][i] * f[j];
        }
    }

    // Обратная подстановка для Rx = qt_f
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = qt_f[i];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= R[i][j] * x[j];
        }
        x[i] /= R[i][i];
    }

    return x;
}
int main() {
    setlocale (LC_ALL, "");
    int n = 1000;
	//матрица А
    vector<vector<double>> a(n, vector<double>(n));
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			a[i][j] = 10 + 0.1 * (i+1) - 0.2 * (j+1);
			a[j][i] = 10 + 0.1 * (j+1) - 0.2 * (i+1);
		}
		a[i][i] = 100;
	}
	//составление f
    vector<double> f(n);
    for (int i = 0; i < n; i++) f[i] = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            f[i] += a[i][j];
        }
    }
	
    // Точное решение (для демонстрации)
    vector<double> x_exact(n, 1.0);

    const int runs = 5;
    vector<double> lu_times, qr_times, lu_errors, qr_errors;

    for (int i = 0; i < runs; ++i) {
        // LU решение
        auto start_lu = high_resolution_clock::now();
        vector<double> x_lu = solveLU(a, f);
        auto end_lu = high_resolution_clock::now();
        lu_times.push_back(duration_cast<microseconds>(end_lu - start_lu).count());
        lu_errors.push_back(calculateError(x_lu, x_exact));

        // QR решение
        auto start_qr = high_resolution_clock::now();
        vector<double> x_qr = solveQR(a, f);
        auto end_qr = high_resolution_clock::now();
        qr_times.push_back(duration_cast<microseconds>(end_qr - start_qr).count());
        qr_errors.push_back(calculateError(x_qr, x_exact));
    }

    // Вычисление средних значений
    double avg_lu_time = accumulate(lu_times.begin(), lu_times.end(), 0.0) / runs;
    double avg_qr_time = accumulate(qr_times.begin(), qr_times.end(), 0.0) / runs;
    double avg_lu_error = accumulate(lu_errors.begin(), lu_errors.end(), 0.0) / runs;
    double avg_qr_error = accumulate(qr_errors.begin(), qr_errors.end(), 0.0) / runs;

    cout << "\nСредние результаты после " << runs << " запусков:" << endl;
    cout << "Метод LU-разложения:" << endl;
    cout << "Среднее время: " << avg_lu_time << " микросекунд" << endl;
    cout << "Средняя погрешность: " << avg_lu_error << endl;

    cout << "\nМетод QR-разложения:" << endl;
    cout << "Среднее время: " << avg_qr_time << " микросекунд" << endl;
    cout << "Средняя погрешность: " << avg_qr_error << endl;

    return 0;

}

