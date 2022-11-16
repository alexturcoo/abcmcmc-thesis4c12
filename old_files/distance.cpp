//This program is a template function for computing the euclidean
//distance between two std::vectors
//taken from https://cplusplus.com/forum/general/209784/
#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

// Computes the distance between two std::vectors
template <typename T>
double vectors_distance(const std::vector<T> &a, const std::vector<T> &b){
    std::vector<double> auxiliary;

    std::transform (a.begin(), a.end(), b.begin(), std::back_inserter(auxiliary), [](T element1, T element2) {return pow((element1-element2),2);});

    return sqrt(std::accumulate(auxiliary.begin(), auxiliary.end(), 0.0));
} // end template vectors_distance

int main() {
    std::vector<double> v1 = {0.1,0.3,0.4};
    std::vector<double> v2 = {-0.1,-0.3, -0.4};

    std::cout << vectors_distance(v1,v2) << "\n";
}
