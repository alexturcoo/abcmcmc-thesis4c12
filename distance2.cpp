#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

//this function I made to normalize vectors
//FEB 23 - THIS IS NOT HOW IM GOING TO NORMALIZE
std::vector<double> normalize_vector(std::vector<double> vec1) {
    std::vector<double> normalized_vec;
    double normalize_value = 0.00;
    int vecSize = vec1.size();
    for (int i = 0; i < vecSize; i++) {
        normalize_value += vec1[i] * vec1[i];
    }
    for (int j = 0; j < vecSize; j++) {
        normalized_vec.push_back(vec1[j]/sqrt(normalize_value));
    }

    // TO TEST THE OUTPUT/PRINT THE VECTOR
    for (int x = 0; x < normalized_vec.size(); x++) {
        std::cout << normalized_vec[x] << ' ' << "\n";
    }


    return normalized_vec;
}

//Feb 23 - New function for normalizing



//this is my function I made to solve euclidean distance
double vectors_distance2(std::vector<double> a, std::vector<double> b) {
    int vecSize = a.size();
    double dist = 0.0;
    for (int i = 0; i < vecSize; i++){
        dist += pow((a[i] - b[i]), 2);
    }

    std::cout << sqrt(dist) << "\n";
    return sqrt(dist);
}

//int main() {
//    std::vector<double> v1 = {0.1,0.3,0.4};
//    std::vector<double> v2 = {-0.1,-0.3, -0.4};

//    std::cout << vectors_distance2(v1,v2) << "\n";

