#include <iostream>
#include <bits/stdc++.h>
#include <vector>
#include <iomanip>
using namespace std;

//Function to find Mean
double Mean(std::vector<double> vec) {
    double sum = 0;
    for (int i=0; i < vec.size(); i++) {
        sum = sum + vec[i];
    }

    return sum/vec.size();
}

//Function to caluclate Standard deviation of given vector
double stdDev(std::vector<double> vec) {
    double sum = 0;
    for (int i = 0; i < vec.size(); i++) {
        sum = sum + (vec[i] - Mean(vec)) * (vec[i] - Mean(vec));
    }

    return sqrt(sum / (vec.size() - 1));
}

//Function to calculate t-test
double tTest(std::vector<double> vec1, double prev_dist) {
    double mean1 = Mean(vec1);
    double sd1 = stdDev(vec1);

    //Equation to get t-statistic for 1 sample
    double t_statistic = abs((mean1 - prev_dist) / (sd1 / sqrt(vec1.size())));

    std::vector<double> areas = {0.005, 0.01, 0.025, 0.05, 0.10};
    std::vector<double> crit_t_values = {3.250, 2.821, 2.262, 1.833, 1.383}; //critical t values 9 degrees of freedom

    if (t_statistic > crit_t_values[0]) {
        std::cout << "P value < 0.005" << "\n";
    }
    else if (t_statistic < crit_t_values[0] && t_statistic > crit_t_values[1]) {
        std::cout << "0.005 < P value < 0.01" << "\n";
    }
    else if (t_statistic < crit_t_values[1] && t_statistic > crit_t_values[2]) {
        std::cout << "0.01 < P value < 0.025" << "\n";
    }
    else if (t_statistic < crit_t_values[2] && t_statistic > crit_t_values[3]) {
        std::cout << "0.025 < P value < 0.05" << "\n";
    }
    else if (t_statistic < crit_t_values[3] && t_statistic > crit_t_values[4]) {
        std::cout << "0.05 < P value < 0.10" << "\n";
    }
    else {
        std::cout << "P value > 0.10" << "\n";
    }

    return abs(t_statistic);


}

//Function to find p value from table above
//Since im doing average of 10, the degrees of freedom would be 10-1
//which is 9 - Using 9 degrees of freedom all the time so crit
//values will be the same all the time
//Larger P value means not statically significant.

//To test the functions
int main() {
    std::vector<double> testvec1 = {10,20,30,40,50};
    double old_dist = 10.00;

    std::cout <<tTest(testvec1, old_dist);

    return 0;
}
