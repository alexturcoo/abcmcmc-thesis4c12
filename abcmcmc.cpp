//#include "myRandomFunctions.c"
#include "functions.cpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <bits/stdc++.h>
#include <fstream> //The fstream library allows us to work with files
#include <vector> // Working with vectors is nicer than arrays, more c++

long int idum = -12345; //FOR RAN1
Ran myran(17); //FOR RAN

//////////////////////////////////////////////
// THE MAIN FUNCTION - WHERE EXECUTION BEGINS/
//////////////////////////////////////////////

int main() {
    
    //srand(time(0)); //Completely random
    srand(10); //Setting seed to test
    //createSeq(100); // Creating the Protein
    og_protein(); // Generating data for SRP40 protein
    std::string mutated_seq = mutateSeqAA(createSeq(100)); // testing first mutation step
    mutateSeqExp(mutated_seq); // testing second mutation step
    return 0;
}

