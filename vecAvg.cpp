#include <iostream>
#include <string>
#include <iterator>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>
#include <numeric>

//this is my function I made to get avg of 10 vectors
std::vector<double> vectors_average(std::vector<vector<double>> a) {
    
    int vecSize = a.size();
    std::vector<double> new_ss_vtr;
    double length = 0.0;
    double num_lcrs = 0.0;
    double avg_entropy = 0.0;

    for (int i = 0; i < vecSize; i++){
        for (int j =0; j < 1; j++){
            //std::cout << avg_entropy << "\n";
            length += a[i][j];
            num_lcrs += a[i][j+1];
            avg_entropy += a[i][j+2];
        }
    }

    length = length/vecSize;
    num_lcrs = num_lcrs/vecSize;
    avg_entropy = avg_entropy/vecSize;

    new_ss_vtr.push_back(length);
    new_ss_vtr.push_back(num_lcrs);
    new_ss_vtr.push_back(avg_entropy);

    //TO TEST THE OUTPUT/PRINT THE VECTOR
    //for (int x = 0; x < new_ss_vtr.size(); x++) {
    //    std::cout << new_ss_vtr[x] << ' ' << "\n";
    //}



    //std::cout << sqrt(dist) << "\n";
    return new_ss_vtr;
}

