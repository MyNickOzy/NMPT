#include "Header.h";

// Функция для создания матрицы A по заданному правилу
MatrixXd createMatrixA(int n) {
    MatrixXd A(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            A(i, j) = 1.0 / (1.0 + 0.2 * (i + 1) + 3.0 * (j + 1));
        }
    }
    return A;
}

// Функция для создания вектора f = A * x*, где x* = (1,1,...,1)^T
VectorXd createVectorF(const MatrixXd& A) {
    int n = A.rows();
    VectorXd x_star = VectorXd::Ones(n);
    return A * x_star;
}

// Метод LU разложения
pair<VectorXd, double> solveLU(const MatrixXd& A, const VectorXd& f) {
    auto start = chrono::high_resolution_clock::now();
    VectorXd x = A.partialPivLu().solve(f);
    auto end = chrono::high_resolution_clock::now();
    return make_pair(x, chrono::duration<double>(end - start).count());
}

// Ручная реализация QR-разложения методом Грама-Шмидта
pair<MatrixXd, MatrixXd> gramSchmidtQR(const MatrixXd& A) {
    int m = A.rows(), n = A.cols();
    MatrixXd Q = MatrixXd::Zero(m, n);
    MatrixXd R = MatrixXd::Zero(n, n);

    for (int j = 0; j < n; ++j) {
        VectorXd v = A.col(j);
        for (int i = 0; i < j; ++i) {
            R(i, j) = Q.col(i).dot(A.col(j));
            v -= R(i, j) * Q.col(i);
        }
        R(j, j) = v.norm();
        Q.col(j) = v / (R(j, j) > 1e-10 ? R(j, j) : 1.0);
    }
    return make_pair(Q, R);
}

// Метод QR разложения (Грамма-Шмидта)
pair<VectorXd, double> solveQR(const MatrixXd& A, const VectorXd& f) {
    auto start = chrono::high_resolution_clock::now();
    auto QR = gramSchmidtQR(A);
    MatrixXd Q = QR.first;
    MatrixXd R = QR.second;
    VectorXd x = R.triangularView<Upper>().solve(Q.transpose() * f);
    auto end = chrono::high_resolution_clock::now();
    return make_pair(x, chrono::duration<double>(end - start).count());
}

// Метод SVD разложения с обнулением малых сингулярных чисел
pair<VectorXd, double> solveSVD(const MatrixXd& A, const VectorXd& f) {
    auto start = chrono::high_resolution_clock::now();

    JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
    VectorXd singular_values = svd.singularValues();
    const double threshold = 1e-12;

    // Обнуляем сингулярные числа меньше threshold
    for (int i = 0; i < singular_values.size(); ++i) {
        if (singular_values(i) < threshold) {
            singular_values(i) = 0.0;
        }
    }

    // Псевдообратная матрица с обнуленными малыми сингулярными числами
    VectorXd x = VectorXd::Zero(A.cols());
    for (int i = 0; i < singular_values.size(); ++i) {
        if (singular_values(i) > 0.0) {
            x += svd.matrixU().col(i).dot(f) / singular_values(i) * svd.matrixV().col(i);
        }
    }

    auto end = chrono::high_resolution_clock::now();
    return make_pair(x, chrono::duration<double>(end - start).count());
}

// Вычисление погрешности ||x* - x|| / ||x*||
double computeError(const VectorXd& x) {
    VectorXd x_star = VectorXd::Ones(x.size());
    return (x_star - x).norm() / x_star.norm();
}

// Вычисление сингулярных чисел и числа обусловленности с обнулением малых значений
tuple<VectorXd, double> computeSVD(const MatrixXd& A) {
    JacobiSVD<MatrixXd> svd(A);
    VectorXd singular_values = svd.singularValues();
    const double threshold = 1e-12;

    // Обнуляем сингулярные числа меньше threshold
    VectorXd modified_singular_values = singular_values;
    for (int i = 0; i < modified_singular_values.size(); ++i) {
        if (modified_singular_values(i) < threshold) {
            modified_singular_values(i) = 0.0;
        }
    }

    // Находим наибольшее и наименьшее ненулевые сингулярные числа
    double d1 = modified_singular_values(0); // Всегда наибольшее
    double dr = d1;
    for (int i = modified_singular_values.size() - 1; i >= 0; --i) {
        if (modified_singular_values(i) > 0.0) {
            dr = modified_singular_values(i);
            break;
        }
    }

    // Если все сингулярные числа нулевые (маловероятно для нашей задачи)
    if (dr == 0.0) dr = threshold;

    double cond = d1 / dr;

    return make_tuple(modified_singular_values, cond);
}

// Функция для вывода всех сингулярных чисел
void printAllSingularValues(const VectorXd& singular_values) {
    cout << "\nAll Singular Values (σ) with values <1e-12 set to 0:\n";
    cout << scientific << setprecision(6);
    for (int i = 0; i < singular_values.size(); ++i) {
        cout << "σ[" << setw(3) << i + 1 << "] = " << singular_values(i);

        // Выделяем наибольшее и наименьшее ненулевые сингулярные числа
        if (i == 0) cout << " (largest)";
        bool is_last_nonzero = true;
        for (int j = i + 1; j < singular_values.size(); ++j) {
            if (singular_values(j) > 0.0) {
                is_last_nonzero = false;
                break;
            }
        }
        if (singular_values(i) > 0.0 && is_last_nonzero) cout << " (smallest non-zero)";
        cout << "\n";
    }
}

int main() {
    const int n = 5; // Размер матрицы
    const int runs = 5; // Количество запусков

    // Создаем матрицу A и вектор f
    MatrixXd A = createMatrixA(n);
    VectorXd f = createVectorF(A);

    // Вычисляем сингулярные числа и число обусловленности
    auto svd_result = computeSVD(A);
    VectorXd singular_values = get<0>(svd_result);
    double cond = get<1>(svd_result);

    // Векторы для хранения результатов
    vector<double> lu_times, qr_times, svd_times;
    vector<double> lu_errors, qr_errors, svd_errors;

    // Проводим измерения
    for (int i = 0; i < runs; ++i) {
        auto lu_result = solveLU(A, f);
        auto qr_result = solveQR(A, f);
        auto svd_result = solveSVD(A, f);

        lu_times.push_back(lu_result.second);
        qr_times.push_back(qr_result.second);
        svd_times.push_back(svd_result.second);

        lu_errors.push_back(computeError(lu_result.first));
        qr_errors.push_back(computeError(qr_result.first));
        svd_errors.push_back(computeError(svd_result.first));
    }

    // Функция для вычисления среднего значения
    auto average = [](const vector<double>& v) {
        return accumulate(v.begin(), v.end(), 0.0) / v.size();
        };

    // Вывод результатов
    cout << scientific << setprecision(6);
    cout << "=== Results for " << n << "x" << n << " matrix ===\n\n";

    cout << "Average Solution Time (s):\n";
    cout << "LU decomposition:  " << average(lu_times) << "\n";
    cout << "QR decomposition:  " << average(qr_times) << "\n";
    cout << "SVD decomposition: " << average(svd_times) << "\n\n";

    cout << "Average Relative Error:\n";
    cout << "LU:  " << average(lu_errors) << "\n";
    cout << "QR:  " << average(qr_errors) << "\n";
    cout << "SVD: " << average(svd_errors) << "\n\n";

    // Подробный вывод информации о числе обусловленности
    cout << "Matrix Condition Number (Cond(A) = σ₁/σᵣ) with <1e-12 set to 0:\n";
    cout << "σ₁ (largest singular value):       " << singular_values(0) << "\n";

    // Находим наименьшее ненулевое сингулярное число
    double smallest_nonzero = 0.0;
    for (int i = singular_values.size() - 1; i >= 0; --i) {
        if (singular_values(i) > 0.0) {
            smallest_nonzero = singular_values(i);
            break;
        }
    }
    cout << "σᵣ (smallest non-zero singular value): " << smallest_nonzero << "\n";
    cout << "Condition number:                  " << cond << "\n";

    // Вывод всех сингулярных чисел
    printAllSingularValues(singular_values);

    return 0;
}