#include "Header.h"
double calculateError(const vector<double>& x) {
    int N = x.size();
    vector<double> x_star(N, 1.0); // ������ ������� x* (������ �� ������)

    // ��������� ������ �������� x - x*
    vector<double> diff(N);
    for (int i = 0; i < N; ++i) {
        diff[i] = x[i] - x_star[i];
    }

    // ��������� �����
    double numerator = euclideanNorm(diff);
    double denominator = euclideanNorm(x_star); // = sqrt(N), �� ��������� ��� ��������

    return numerator / denominator;
}
