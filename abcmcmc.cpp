//#include "myRandomFunctions.c"
#include "functions.cpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <bits/stdc++.h>
#include <fstream> //The fstream library allows us to work with files
#include <vector> // Working with vectors is nicer than arrays, more c++

// THIS IS A TEST TO SEE IF I CAN PUSH AND PULL FROM GITHUB - WORKS
// PUT THINGS YOURE GONNA DEFINE HERE
// 21 amino acids because of 0 at the end which is not an amino acid
#define numAA 20

long int idum = -12345; //FOR RAN1
Ran myran(17); //FOR RAN


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
    
    //std::cout << protein << "\n" << "\n" ;
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
    //std::cout << srp40 << "\n" << "\n" << length << "\n" << "\n"<< num_lcr << "\n" << "\n" ;
    return 0;
}

/*THIS FUNCTION WILL CHOOSE A RANDOM POISSON DEVIATE WITH A MEAN
EQUAL TO (MUTATION RATE * # OF SITES) - WHATEVER THE VALUE OF THE
RANDOM DEVIATE, CHOOSE THAT MANY MUTATIONS IN THE SEQUENCE AT
RANDOM*/

std::string mutateSeqAA(std::string simulated_protein){

    float mutation_rate = 0.14 ; // just setting mutation rate here now
    float avg_mutations = mutation_rate * simulated_protein.length() ;
    Expondev myexpdev(avg_mutations, 100); //creating an instance
    float poisson_deviate_num_mutations = myexpdev.dev() ; //Choose random poisson deviate with mean = avg_mutations
    std::cout << simulated_protein << "\n" << poisson_deviate_num_mutations << "\n" ; 
    //This for loop will give 5 random sites mutations
    for (int i = 0; i < poisson_deviate_num_mutations; i++) {

        char aminoAcids[numAA] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 
                                   'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };
        
        char random_AA = aminoAcids[rand() % numAA]; // sets up the random amino acid, same used in first function to createSeq
        int random_site = rand() % simulated_protein.length(); // Gives random number in range (0 - length of protein sequence)
        simulated_protein[random_site] = random_AA; // indexes the simulated protein at a random spot and replaces the existing AA with a new random one
    }

    std::cout << "after mutateSeqAA: " << simulated_protein << "\n" ;
    return simulated_protein ;
}

/*THIS FUNCTION WILL MUTATE THE SIMULATED PROTEIN SEQUENCE
BY CHOOSING A RANDOM EXPONENTIAL DEVIATE (WITH MEAN = MUTATION RATE)
FOR EACH AMINO ACID IN THE SEQUENCE AND SUBSEQUENTLY SELECTING
THE AMINO ACID WITH THE LOWEST NUMBER (QUICKEST TO MUTATE) AND
MUTATING IT RANDOMLY, THIS IS DONE SUCCESSIVELY TO PRODUCE A 
PROTEIN AND CREATE A VECTOR OF VALUES SIMILAR TO ABOVE (SEP 21)     
SEP 26 - PROCESS CHANGE, THIS FUNCTION WILL BE USED FOR AMINO ACID EXPANSION */

std::string mutateSeqExp(std::string simulated_protein){

    std::vector<double> exp_deviates_vtr ; // Creating a vector to hold the values of the deviates

    // Traversing the string
    for (int i = 0; i < simulated_protein.length(); i++) {

        int counter = 1 ;

        //Code to scan back and forth to find repeats
        if (simulated_protein[i] != simulated_protein[i+1] && simulated_protein[i] != simulated_protein[i-1]) {
            float mutation_rate = 0.14 ; //mutation rate == alpha
            float beta = mutation_rate*1 ; // 1 will always be used here because the length if no repeats is 1
            double deviate = myran.doub();//here we choose exp_deviate(mean of beta)
            exp_deviates_vtr.push_back(deviate) ; //Here we are storing the exponential deviates
        } else {
            int x = 1 ;
            int y = 1 ;
            
            //Be careful in these while loops, for i-y, when i is 0
            //and y is 1, how does it not throw error
            while (simulated_protein[i] == simulated_protein[i + x]) {
                counter += 1 ;
                x++;
            } 

            while (simulated_protein[i] == simulated_protein[i - y]) {
                counter += 1 ;
                y++;
            }

            float mutation_rate = 0.14 ;
            float beta = mutation_rate * counter ;
            int deviate = rand() % 10;
            exp_deviates_vtr.push_back(beta);
        }
    }
    // THIS IS JUST TO PRINT THE VECTOR
    for (int x = 0; x < exp_deviates_vtr.size(); x++) {
        std::cout << exp_deviates_vtr[x] << ' ';
    }
    
    std::cout << "after mutateSeqEXP: " << simulated_protein << "\n" << "\n" ;
    return simulated_protein;
}

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

