#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <bits/stdc++.h>
#include <fstream> //The fstream library allows us to work with files
#include "myRandomFunctions.c"
#include <vector> // Working with vectors is nicer than arrays, more c++

// THIS IS A TEST TO SEE IF I CAN PUSH AND PULL FROM GITHUB - WORKS
// PUT THINGS YOURE GONNA DEFINE HERE
// 21 amino acids because of 0 at the end which is not an amino acid
#define numAA 20

////////////////////////////////////////////////////////////////
// FIRST FUNCTION TAKES AN INTEGER VALUE AND GENERATES A RANDOM/
// AMINO ACID SEQUENCE OF THAT LENGTH                          / 
////////////////////////////////////////////////////////////////

std::string createSeq(int n){

    char aminoAcids[numAA] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 
                                   'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };

    std::string protein = "";
    for (int i = 0; i < n; i++){
        protein = protein + aminoAcids[rand() % numAA];}
    
    // std::cout << protein << "\n" << "\n" ;
    // printf("%s", protein.c_str());
    return protein;
}

////////////////////////////////////////////////////////////////////
// SECOND FUNCTION READS IN SRP40 PROTEIN AND CREATES A VECTOR WITH/
// CHARACTERISTICS - THIS IS EITHER S or D - THE ORIGINAL DATA     /
////////////////////////////////////////////////////////////////////

int og_protein(){

    std::string srp40 = "MASKKIKVDEVPKLSVKEKEIEEKSSSSSSSSSSSSSSSSSSSSSSSSSGESSSSSSSSSSSSSSDSSDS\
SDSESSSSSSSSSSSSSSSSDSESSSESDSSSSGSSSSSSSSSDESSSESESEDETKKRARESDNEDAKE\
TKKAKTEPESSSSSESSSSGSSSSSESESGSESDSDSSSSSSSSSDSESDSESDSQSSSSSSSSDSSSDS\
DSSSSDSSSDSDSSSSSSSSSSDSDSDSDSSSDSDSSGSSDSSSSSDSSSDESTSSDSSDSDSDSDSGSS\
SELETKEATADESKAEETPASSNESTPSASSSSSANKLNIPAGTDEIKEGQRKHFSRVDRSKINFEAWEL\
TDNTYKGAAGTWGEKANEKLGRVRGKDFTKNKNKMKRGSYRGGSITLESGSYKFQD";

    //creating fasta file of protein string to use in segA
    std::ofstream MyFile("srp40.fasta"); //Create and open a fasta file
    MyFile << ">srp40_saccharomyces_cerevisiae\n" << srp40; //Write to the file
    MyFile.close(); //Close the file

    //finding the LENGTH + #ofLCRs
    int length = srp40.length();
    system("segA srp40.fasta 15 1.9 2.2 -l | grep '>' | wc -l >> numlcr.txt"); //Write the results of counting the LCRs to a txt file
    std::ifstream myreadfile("numlcr.txt"); //Declaring fstream variable and connecting it to a stream object by opening the file
    int num_lcr; //Declare int variable which holds the contents of the 'numlcr.txt' file
    myreadfile >> num_lcr; //Pipe file's content into stream

    //finding the average entropy of the LCRs
    //int avg_entropy_lcr =


    //TO TEST THE OUTPUT
    std::cout << srp40 << "\n" << "\n" << length << "\n" << "\n"<< num_lcr << "\n" << "\n" ;
    return 0;
}

////////////////////////////////////////////////////////////////////////
// THIS FUNCTION WILL MUTATE THE SIMULATED PROTEIN SEQUENCE
// BY CHOOSING A RANDOM EXPONENTIAL DEVIATE (WITH MEAN = MUTATION RATE)
// FOR EACH AMINO ACID IN THE SEQUENCE AND SUBSEQUENTLY SELECTING
// THE AMINO ACID WITH THE LOWEST NUMBER (QUICKEST TO MUTATE) AND
// MUTATING IT RANDOMLY, THIS IS DONE SUCCESSIVELY TO PRODUCE A 
// PROTEIN AND CREATE A VECTOR OF VALUES SIMILAR TO ABOVE (SEP 21)     /
////////////////////////////////////////////////////////////////////////

std::string mutateSeq(std::string simulated_protein){

    std::vector<double> exp_deviates_vtr ; // Creating a vector to hold the values of the deviates

    // Traversing the string
    for (int i = 0; i < simulated_protein.length(); i++) {

        // Code to apply random exponential deviates to AA
        int deviate =  rand();
        exp_deviates_vtr.push_back(deviate);
    }

    // THIS IS JUST TO PRINT THE VECTOR
    for (int x = 0; x < exp_deviates_vtr.size(); x++) {
        std::cout << exp_deviates_vtr[x] << ' ';
    }
    

    return simulated_protein;
}

//////////////////////////////////////////////
// THE MAIN FUNCTION - WHERE EXECUTION BEGINS/
//////////////////////////////////////////////

int main() {
    
    srand(18);
    //createSeq(100); // Creating the Protein
    og_protein(); // Generating data for SRP40 protein
    mutateSeq(createSeq(10)); // Mutating the simulated Sequence
    return 0;
}

