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
std::vector<double> normalize_observed(std::vector<double> vec1) {
    std::vector<double> normalized_vec;
    int vecSize = vec1.size();
    int length = vec1[0];
    int num_lcrs = vec1[1];
    double avg_entropy = vec1[2];

    //setting the max and mins
    int max_length = length*1.5;
    int min_length = length*0.5;
    double max_lcrs = length/8; //8 is the hypothetical length of the shortest LCR, this will give you a maximal total number of LCRs
    int min_lcrs = 0;
    double max_entropy = 4.3;

    //pushing back the normalized values
    normalized_vec.push_back((length - (max_length-length) + min_length)/406);
    normalized_vec.push_back(num_lcrs/max_lcrs);
    normalized_vec.push_back(avg_entropy/max_entropy);

    // TO TEST THE OUTPUT/PRINT THE VECTOR
    for (int x = 0; x < normalized_vec.size(); x++) {
        std::cout << normalized_vec[x] << ' ' << "\n";
    }

    return normalized_vec;

}

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

int main() {
    std::vector<double> v1 = {406,3,1.47667};
    normalize_observed(v1);
}

