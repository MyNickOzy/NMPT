#include "Header.h"
// Функция для форматированного вывода числа с высокой точностью
string formatNumber(double value, int precision = 15) {
    ostringstream oss;
    oss << scientific << setprecision(precision) << value;
    string s = oss.str();

    // Убираем лишние нули в экспоненциальной записи
    size_t e_pos = s.find('e');
    if (e_pos != string::npos) {
        size_t dot_pos = s.find('.');
        if (dot_pos != string::npos) {
            // Удаляем лишние нули после точки
            size_t last_non_zero = s.find_last_not_of('0', e_pos - 1);
            if (last_non_zero != string::npos && last_non_zero > dot_pos) {
                s.erase(last_non_zero + 1, e_pos - last_non_zero - 1);
            }
        }
    }

    // Если экспонента равна +00, заменяем на обычную запись
    if (s.find("e+00") != string::npos) {
        s = s.substr(0, s.find("e+00"));
        // Удаляем точку, если после неё нет цифр
        if (s.find('.') == s.length() - 1) {
            s = s.substr(0, s.length() - 1);
        }
    }

    return s;
}
int main() {
    const int N = 1000;  // Размер матрицы
    const int runs = 5;  // Количество запусков для усреднения

    // Создаем матрицу A и вектор f (однократно)
    vector<vector<double>> A = createMatrixA(N);
    vector<double> f = createVectorF(A);

    // Переменные для накопления результатов
    double total_time_LU = 0.0;
    double total_error_LU = 0.0;
    double total_time_QR = 0.0;
    double total_error_QR = 0.0;

    cout << "Running " << runs << " iterations for N = " << N << endl;

    for (int run = 0; run < runs; ++run) {
        // Решение методом LU
        auto start_LU = high_resolution_clock::now();
        vector<double> x_LU = solveWithLU(A, f);
        auto stop_LU = high_resolution_clock::now();
        auto duration_LU = duration_cast<nanoseconds>(stop_LU - start_LU);

        // Решение методом QR
        auto start_QR = high_resolution_clock::now();
        vector<double> x_QR = solveWithQR(A, f);
        auto stop_QR = high_resolution_clock::now();
        auto duration_QR = duration_cast<nanoseconds>(stop_QR - start_QR);

        // Вычисление погрешностей
        double error_LU = calculateError(x_LU);
        double error_QR = calculateError(x_QR);

        // Накопление результатов
        total_time_LU += duration_LU.count();
        total_error_LU += error_LU;
        total_time_QR += duration_QR.count();
        total_error_QR += error_QR;
    }

    // Вычисление средних значений
    double avg_time_LU = total_time_LU / runs;
    double avg_error_LU = total_error_LU / runs;
    double avg_time_QR = total_time_QR / runs;
    double avg_error_QR = total_error_QR / runs;

    // Вывод результатов с высокой точностью
    cout << "\nAverage results after " << runs << " runs:" << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "Method\t\tTime (ns)\t\tRelative Error (L)" << endl;
    cout << "------------------------------------------------------------" << endl;
    cout << "LU\t\t" << formatNumber(avg_time_LU) << "\t" << formatNumber(avg_error_LU) << endl;
    cout << "QR\t\t" << formatNumber(avg_time_QR) << "\t" << formatNumber(avg_error_QR) << endl;
    cout << "------------------------------------------------------------" << endl;

    return 0;
}
