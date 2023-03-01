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
/*std::vector<double> normalize_vector(std::vector<double> vec1) {
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
}*/

//Feb 23 - New function for normalizing
std::vector<double> normalize_vector(std::vector<double> vec1) {
    std::vector<double> normalized_vec;
    int srp40_length = 406;
    int vecSize = vec1.size();
    float length = vec1[0];
    float num_lcrs = vec1[1];
    double avg_entropy = vec1[2];

    //setting the max and mins
    float max_length = srp40_length*1.5;
    float min_length = srp40_length*0.5;

    double max_lcrs = srp40_length/8; //8 is the hypothetical length of the shortest LCR, this will give you a maximal total number of LCRs
    double max_entropy = 4.3;

    //pushing back the normalized values
    normalized_vec.push_back((length - srp40_length + min_length)/srp40_length);
    normalized_vec.push_back(num_lcrs/max_lcrs);
    normalized_vec.push_back(avg_entropy/max_entropy);

    // TO TEST THE OUTPUT/PRINT THE VECTOR
    /*for (int x = 0; x < normalized_vec.size(); x++) {
        std::cout << normalized_vec[x] << ' ' << "\n";
    }*/

    return normalized_vec;

}

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

//FOR TESTING THIS FILE
/*int main() {
    std::vector<double> v1 = {500,6,2.678};
    normalize_vector(v1);
}*/

