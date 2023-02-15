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

Ran myran(time(NULL)); //We will use 21 as the random seed right now, used in Poissondev too

////////////////////////////////////////////////////////////////
// FIRST FUNCTION TAKES AN INTEGER VALUE AND GENERATES A RANDOM/
// AMINO ACID SEQUENCE OF THAT LENGTH                          / 
////////////////////////////////////////////////////////////////
std::string createSeq(int n){

    char aminoAcids[numAA] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 
                                   'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };

    std::string protein = "";
    for (int i = 0; i < n; i++){
        protein += aminoAcids[myran.int64() % numAA];} //this rand() % 20 means in the range 0-19
    
    //std::cout << protein << "\n" << "\n" ;
    return protein;
}

///////////////////////////////////////////////////
/////// FUNCTION TO GENERATE NORMAL DEVIATES //////
///////////////////////////////////////////////////
double getNormalDev(double mu, double stdev) {
    Normaldev mynorm(mu, stdev, myran.int64());
    double dev = mynorm.dev();
    //std::cout << dev << "\n";
    return dev;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////// Changed here down - Brian - Feb 12 Xander made Changes, trying to understand the new stuff ////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::string mutateSeqExpBG(std::string simulated_protein, float mutation_rate, float indel_rate){ //took out iterates as a parameter

    //Setting up the vectors
    std::vector<double> mut_dev;
    std::vector<double> ind_dev;

    char aminoAcids[20] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E',
        'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };
                                      
    std::cout << "before mutateseqEXP:\t" << simulated_protein << "\n"; // Initially printing the non-mutated strin.
    int len = simulated_protein.length();

    //ASSIGN DEVIATES FOR MUTATION
    for (int i = 0; i < len; i++) {
        float beta1 = mutation_rate ; 
        Expondev myexp(beta1,myran.int64());
        double deviate = myexp.dev(); // choose exp_deviate(mean of beta)
        mut_dev.push_back(deviate);
    }
    // Traverse the string and generate deviates FOR INDELS - this seems to be
    // working, counter is showing correct values for repeats -Feb 15
    for (int i = 0; i < len; i++) {
        int counter = 1 ;
        //Code to scan back and forth to find repeats - this i-1 will give -1 on the first iteration though - check this in while loop too for scanning backwards
        if (simulated_protein[i] != simulated_protein[i+1] && simulated_protein[i] != simulated_protein[i-1]) {
            float beta2 = indel_rate ; // the length if no repeats is 1
            Expondev myexp(beta2,myran.int64());
            double deviate = myexp.dev();
            ind_dev.push_back(deviate);  
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
            float beta3 = indel_rate * counter ; // trying to see if the # of repeats plays a role, multiply by counter
            Expondev myexp(beta3,myran.int64());
            double deviate = myexp.dev();
            ind_dev.push_back(deviate);
        }
        //std::cout << counter << "\n"; //Checking if counter is finding repeats
    }

    // Working with the vector of structs, why we doing 0.5*len
    // FEB 15 - WORKING WITH VECTORS NOW 
    for(int iter=0; iter<0.5*len; iter++) { // throw down 0.5 mut/site
        //selecting the lowest deviate from both vectors
        double minExpDev=mut_dev[0];
        int minPosition=0;
        int minType=0; // type=0 for mut, 1 for ins, 2 for del
        for (int i = 0; i < len; i++) {
            if(minExpDev > mut_dev[i]) {
                minExpDev=mut_dev[i];
                minPosition=i; minType=0;
            }
            if(minExpDev > ind_dev[i]) {
                minExpDev=ind_dev[i];
                minPosition=i; minType=1;
                if(myran.doub() < 0.5) { minType=2; }
            }
        }
        std::cout << "minpos" << minPosition << "\n";
        std::cout << "mindev" << minExpDev << "\n";
        std::cout << "mintype" << minType << "\n";
        std::cout << "iter" << iter << "\n";

        //This is where base changes, insertions, deletions occur
        if(minType==0) { // put in a base change
            simulated_protein[minPosition]=aminoAcids[myran.int64() % numAA];
            float beta1 = mutation_rate ; 
            Expondev myexp(beta1,myran.int64());
            double deviate = myexp.dev();
            mut_dev[minPosition]=deviate;//swap deviate for base change
            //Printing to see mutated protein afer base change//
            std::cout << "\n" << "after basechange:\t" << simulated_protein << "\n" << "\n" ;
            // THIS IS JUST TO PRINT THE VECTOR OF DEVIATES FOR MUTATION TO SEE IF IT CHANGED
            for (int x = 0; x < mut_dev.size(); x++) {
                std::cout << mut_dev[x] << ' ';
            }
            for (int x = 0; x < ind_dev.size(); x++) {
                std::cout << ind_dev[x] << ' ';
            }
        } else { // put in an indel
            if(minType==1) { // use insertion
                if(len < 1.5*len) {
                    len++;
                    simulated_protein.insert(minPosition+1, 1, simulated_protein[minPosition]); // insert repeat beside minPosition
                    float beta1 = mutation_rate ; 
                    Expondev myexp(beta1,myran.int64());
                    double deviate = myexp.dev();
                    mut_dev.insert(mut_dev.begin() + minPosition+1, deviate); //Adding a new deviate for the insertion
                    //Printing to see mutated protein after insertion//
                    std::cout << "\n" << "after insertion:\t" << simulated_protein << "\n" << "\n" ;
                    // THIS IS JUST TO PRINT THE VECTOR OF DEVIATES FOR MUTATION TO SEE IF IT CHANGED
                    for (int x = 0; x < mut_dev.size(); x++) {
                        std::cout << mut_dev[x] << ' ';
                    }
                    for (int x = 0; x < ind_dev.size(); x++) {
                        std::cout << ind_dev[x] << ' ';
                    }
                }
            } else { // use deletion
                if(len > 0.5*len) {
                    len--;
                    simulated_protein.erase(minPosition, 1); //delete amino acid at minPosition
                    mut_dev.erase(mut_dev.begin() + minPosition); //Deleting deviate for the deletion
                    //Printing to see mutated protein after deletion//
                    std::cout << "\n" << "after deletion:\t" << simulated_protein << "\n" << "\n" ;
                    // THIS IS JUST TO PRINT THE VECTOR OF DEVIATES FOR MUTATION TO SEE IF IT CHANGED
                    for (int x = 0; x < mut_dev.size(); x++) {
                        std::cout << mut_dev[x] << ' ';
                    }
                    for (int x = 0; x < ind_dev.size(); x++) {
                        std::cout << ind_dev[x] << ' ';
                    }
                }
            }
        }
        
        // check what is around the change
        // need to look down one, the site and up one
        // So here were checking if the insertion,
        // or deletion affected the landscape of the DNA sequence.
        // Did it create a repeat, inturrupt a repeat, etc.
        /*if (simulated_protein[minPosition] == simulated_protein[minPosition+1] || simulated_protein[minPosition] == simulated_protein[minPosition-1]) { // part of a repeat
            int counter = 1;
            int x = 1;
            int y = 1;
            //Look backward
            while (simulated_protein[minPosition] == simulated_protein[minPosition-y]) {
                counter += 1;
                y++;
            }
            int start = minPosition - y; //Getting start position of repeats
            std::cout << "start" << start << "\n";
            //count forward
            while (simulated_protein[minPosition] == simulated_protein[minPosition + x]) { 
                counter += 1; 
                x++; 
            }
            int end = minPosition + x; //Getting end position of repeats
            std::cout << "end" << end << "\n";
            float beta3 = indel_rate * counter ;
            for(int j=start; j<=end; j++) {
                Expondev myexp(beta3,myran.int64());
                double deviate = myexp.dev();
                ind_dev[j]=deviate;
            }
        } else { // not part of a repeat
            float beta3 = indel_rate;
            Expondev myexp(beta3,myran.int64());
            double deviate = myexp.dev();
            ind_dev[minPosition]=deviate;
        }*/
    }

    std::cout << "\n" << "after mutateSeqEXP:\t" << simulated_protein << "\n" << "\n" ;
    return simulated_protein; 
}

///// TRYING TO DEBUG THIS FILE WITH THIS MAIN FUNCTION /////
int main() {
    for (int j = 0; j < 1; j++){
        std::string simulated_protein = createSeq(10);
        double mutation_rate = 0.14;
        double indel_rate = 0.14;
        std::string mutated_protein = mutateSeqExpBG(simulated_protein, mutation_rate, indel_rate);
        simulated_protein = mutated_protein;
    }
}
