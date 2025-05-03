#include "Header.h"
double calculateError(const vector<double>& x) {
    int N = x.size();
    vector<double> x_star(N, 1.0); // Точное решение x* (вектор из единиц)

    // Вычисляем вектор разности x - x*
    vector<double> diff(N);
    for (int i = 0; i < N; ++i) {
        diff[i] = x[i] - x_star[i];
    }

    // Вычисляем нормы
    double numerator = euclideanNorm(diff);
    double denominator = euclideanNorm(x_star); // = sqrt(N), но вычисляем для общности

    return numerator / denominator;
}
