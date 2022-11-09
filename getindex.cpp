//This function was taken from a website
//geeksforkeeks.org/how-to-find-index-of-a-given-element etc

//c++ program to find the index
//of an element in a vector

#include <bits/stdc++.h>
using namespace std;

//Function to print the index of an element
int getIndex(vector<double> v, double k){

    auto it = find(v.begin(), v.end(), k);

    // If element was found
    if (it != v.end()){

        //calculating the index of k
        int index = it - v.begin();
        return index;
    }
    else {
        //if the element is not present in the vector
        return -1;
    }
}
