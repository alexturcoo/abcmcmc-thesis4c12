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

//Function to get the equation of a line and estimate the
//probability based on the t-statistic
double eqn_line(double y2,double y1,double x2,double x1,double tstat) {
    double m = (y2 - y1) / (x2 - x1);
    double b = y2 - (m*x2);
    double x = (tstat - b) / m;

    return x;
}

//ReWriting function to obtain estimate of pvalue from linear
//interpolation
double tTest(std::vector<double> vec1, double prev_dist) {
    double mean1 = Mean(vec1);
    double sd1 = stdDev(vec1);

    //Equation to get t-statistic for 1 sample
    double t_statistic = abs((mean1 - prev_dist) / (sd1 / sqrt(vec1.size())));

    std::vector<double> areas = {0.5, 0.25, 0.2, 0.15, 0.10, 0.05, 0.025, 0.01, 0.005, 0.001, 0.0005}; //this is like y
    std::vector<double> crit_t_values = {0.000, 0.703, 0.883, 1.100, 1.383, 1.833, 2.262, 2.821, 3.250, 4.297, 4.781}; //critical t values 9 degrees of freedom - this is like x

    int minPosition=0;
    for (int i = 0; i < areas.size(); i++) {
        if (t_statistic > crit_t_values[i]) {
            minPosition = i;
        }
    }
    int maxPosition = minPosition + 1;

    double pval_est = eqn_line(crit_t_values[minPosition], crit_t_values[maxPosition], areas[minPosition], areas[maxPosition], t_statistic);

    return pval_est;
    return minPosition; 
    return abs(t_statistic);
}


//Function to find p value from table above
//Since im doing average of 10, the degrees of freedom would be 10-1
//which is 9 - Using 9 degrees of freedom all the time so crit
//values will be the same all the time
//Larger P value means not statically significant.

//To test the functions
int main() {
    std::vector<double> testvec1 = {10,20,30,40,50,60,70,80,90,100};
    double old_dist = 27.00;

    std::cout << tTest(testvec1, old_dist);

    return 0;
}
