#include "functions.cpp"
#include "getindex.cpp"
#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <string>
#include <ctime>
#include <algorithm> //this is to get min element stuff, cool library
#define numAA 20
using namespace std;


Ran myran(21); //We will use 21 as the random seed right now, used in Poissondev too
/*THIS FUNCTION WILL CHOOSE A RANDOM POISSON DEVIATE WITH A MEAN
EQUAL TO (MUTATION RATE * # OF SITES) - WHATEVER THE VALUE OF THE
RANDOM DEVIATE, CHOOSE THAT MANY MUTATIONS IN THE SEQUENCE AT
RANDOM*/

std::string mutateSeqAA(std::string simulated_protein){

    double mutation_rate = 0.14 ; // just setting mutation rate here now
    double avg_mutations = mutation_rate * simulated_protein.length() ;
    Poissondev mypoiss(avg_mutations, 21); //creating an instance
    int poisson_deviate_num_mutations = mypoiss.dev() ; //Choose random poisson deviate with mean = avg_mutations
    //std::cout << poisson_deviate_num_mutations << "\n"
    std::cout << "before mutateseqAA:\t" << simulated_protein << "\n" ; 
    //This for loop will give 5 random sites mutations
    for (int i = 0; i < poisson_deviate_num_mutations; i++) {

        char aminoAcids[numAA] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 
                                   'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };
        
        char random_AA = aminoAcids[myran.int64() % numAA]; // sets up the random amino acid, same used in first function to createSeq
        int random_site = myran.int64() % simulated_protein.length(); // Gives random number in range (0 - length of protein sequence)
        simulated_protein[random_site] = random_AA; // indexes the simulated protein at a random spot and replaces the existing AA with a new random one
    }

    std::cout << "after mutateSeqAA:\t" << simulated_protein << "\n" ;
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
            float beta = mutation_rate ; // 1 will always be used here because the length if no repeats is 1
            Expondev myexp(beta,myran.int64());
            double deviate = myexp.dev();//here we choose exp_deviate(mean of beta)
            exp_deviates_vtr.push_back(deviate) ; //Here we are storing the exponential deviates
        } else {
            int x = 1 ;
            int y = 1 ;
            
            //Be careful in these while loops, for i-y, when i is 0
            //and y is 1, how does it not throw error
            //Looking forward for repeats
            while (simulated_protein[i] == simulated_protein[i + x]) {
                counter += 1 ;
                x++;
            } 
            //Looking backwards for repeats
            while (simulated_protein[i] == simulated_protein[i - y]) {
                counter += 1 ;
                y++;
            }

            float mutation_rate = 0.14 ;
            float beta = mutation_rate * counter ;
            Expondev myexp(beta,myran.int64());
            double deviate = myexp.dev();
            exp_deviates_vtr.push_back(deviate);
        }
    }

    //selecting the lowest deviate and inserting or deleting am acid
    double min = *min_element(exp_deviates_vtr.begin(), exp_deviates_vtr.end());
    int position = getIndex(exp_deviates_vtr, min);
    char aa_index = simulated_protein[position];
    float random_number = myran.doub();

    //Inserting or deleting a repeat
    if (random_number < 0.5){
        simulated_protein.erase(position, 1);
    } else {
        simulated_protein.insert(position+1,1, aa_index);
    }

    //printing out this stuff to check its working
    //std::cout << min << "\n" << position  << "\n" << aa_index << "\n" << random_number << "\n";

    // THIS IS JUST TO PRINT THE VECTOR
    //for (int x = 0; x < exp_deviates_vtr.size(); x++) {
    //    std::cout << exp_deviates_vtr[x] << ' ';
    //}
    
    std::cout << "\n" << "after mutateSeqEXP:\t" << simulated_protein << "\n" << "\n" ;
    return simulated_protein;
}

int main() {
    std::string x = mutateSeqAA("MKNHCHKISAKHHHHAM");
    mutateSeqExp(x); 
}

