#include <vector> 
#include "CG2.h"

using namespace std;
using namespace PANSFEM2;

int main() {
    vector<double> A = { 1, 2, 3, 4 };
    vector<double> b = { 5, 3 };
    vector<double> x(2, 0.0);

    ScalingBiCGSTAB2(A.data(), b.data(), x.data(), 2, 1000, 1.0e-5, true);

    for (auto xi : x) {
        std::cout << xi << std::endl;
    }
}