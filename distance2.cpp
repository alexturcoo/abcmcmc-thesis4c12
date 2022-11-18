#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

//this is my function I made to solve euclidean distance
double vectors_distance2(std::vector<double> a, std::vector<double> b) {
    int vecSize = a.size();
    double dist = 0.0;
    for (int i = 0; i < vecSize; i++){
        dist += pow((a[i] - b[i]), 2);
    }

    //std::cout << sqrt(dist) << "\n";
    return sqrt(dist);
}

//int main() {
//    std::vector<double> v1 = {0.1,0.3,0.4};
//    std::vector<double> v2 = {-0.1,-0.3, -0.4};

//    std::cout << vectors_distance2(v1,v2) << "\n";
//}
