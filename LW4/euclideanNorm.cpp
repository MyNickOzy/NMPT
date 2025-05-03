#include "Header.h"
double euclideanNorm(const vector<double>& v) {
    double norm = 0.0;
    for (double val : v) {
        norm += val * val;
    }
    return sqrt(norm);
}