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


/*THIS FUNCTION WILL MUTATE THE SIMULATED PROTEIN SEQUENCE
BY CHOOSING A RANDOM EXPONENTIAL DEVIATE (WITH MEAN = MUTATION RATE)
FOR EACH AMINO ACID IN THE SEQUENCE AND SUBSEQUENTLY SELECTING
THE AMINO ACID WITH THE LOWEST NUMBER (QUICKEST TO MUTATE) AND
MUTATING IT RANDOMLY, THIS IS DONE SUCCESSIVELY TO PRODUCE A 
PROTEIN AND CREATE A VECTOR OF VALUES SIMILAR TO ABOVE (SEP 21)     
SEP 26 - PROCESS CHANGE, THIS FUNCTION WILL BE USED FOR AMINO ACID EXPANSION
UPDATE NOV 15 - WE ARE MAKING THIS NEW FUNCTION THAT ASSIGNS DEVIATES
BASED ON BOTH MUTATION RATES AND INDEL RATES AND THEN OUT OF BOTH VECTORS
WE WILL FIND THE LOWEST AND EITHER MUTATE IT OR DO AN INS/DEL DEPENDING
ON WHICH VECTOR THE LOWEST DEVIATE CAME FROM*/

std::string mutateSeqExp(std::string simulated_protein, float mutation_rate, float indel_rate){

    // Setting up the vectors
    std::vector<double> exp_deviates_vtr_ind;  // Creating a vector to hold the values of the deviates for indel rate
    std::vector<double> exp_deviates_vtr_mut; // Creating a vector to hold the values of the deviates for mutation rate
    std::vector<double> smallest_vtr; // Creating a vector to store the smallest element of each of the 2 vectors
                                      
    std::cout << "before mutateseqEXP:\t" << simulated_protein << "\n"; // Initially printing the non-mutated strin.

    // Mutation and indel rate set here now - dont need to set here
    //float mutation_rate = 0.14;
    //float indel_rate = 0.14;

    // First loop will assign deviates based on mutation rates
    for (int i = 0; i < simulated_protein.length(); i++) {

        float beta1 = mutation_rate ; // 1 will always be used here because the length if no repeats is 1
        Expondev myexp(beta1,myran.int64());
        double deviate = myexp.dev();//here we choose exp_deviate(mean of beta)
        exp_deviates_vtr_mut.push_back(deviate) ; //Here we are storing the exponential deviates
    }
    // Traversing the string
    for (int i = 0; i < simulated_protein.length(); i++) {

        int counter = 1 ;

        //Code to scan back and forth to find repeats
        if (simulated_protein[i] != simulated_protein[i+1] && simulated_protein[i] != simulated_protein[i-1]) {
            float beta2 = indel_rate ; // 1 will always be used here because the length if no repeats is 1
            Expondev myexp(beta2,myran.int64());
            double deviate = myexp.dev();//here we choose exp_deviate(mean of beta)
            exp_deviates_vtr_ind.push_back(deviate) ; //Here we are storing the exponential deviates
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

            float beta3 = indel_rate * counter ;
            Expondev myexp(beta3,myran.int64());
            double deviate = myexp.dev();
            exp_deviates_vtr_ind.push_back(deviate);
        }
    }

    //selecting the lowest deviate from both vectors
    double min_mut = *min_element(exp_deviates_vtr_mut.begin(), exp_deviates_vtr_mut.end());
    double min_ind = *min_element(exp_deviates_vtr_ind.begin(), exp_deviates_vtr_ind.end());

    // Append the 2 minimums to the new vector
    smallest_vtr.push_back(min_mut);
    smallest_vtr.push_back(min_ind);

    // Get the smallest of the small numbers
    double smallest_num = *min_element(smallest_vtr.begin(), smallest_vtr.end());

    // Getting the index of the smallest of the small numbers
    int position = getIndex(smallest_vtr, smallest_num);

    // If index = 0 mutation, If index = 1, indel
    if (position == 0){

        for (int i = 0; i < 1; i++) {

            char aminoAcids[20] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 
                                   'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };
        
            char random_AA = aminoAcids[myran.int64() % numAA]; // sets up the random amino acid, same used in first function to createSeq
            int position2 = getIndex(exp_deviates_vtr_mut, min_mut);
            simulated_protein[position2] = random_AA; // indexes the simulated protein at a random spot and replaces the existing AA with a new random one
            //std::cout << position2 << "\n";
        }

        // THIS IS JUST TO PRINT THE VECTOR OF DEVIATES FOR MUTATION
            for (int x = 0; x < exp_deviates_vtr_mut.size(); x++) {
                std::cout << exp_deviates_vtr_mut[x] << ' ';
            }

    } else {
        
        int position3 = getIndex(exp_deviates_vtr_ind, min_ind);
        char aa_index = simulated_protein[position3];
        float random_number = myran.doub();

        //Inserting or deleting a repeat
        if (random_number < 0.5){
        simulated_protein.erase(position3, 1);
        } else {
        simulated_protein.insert(position3+1,1, aa_index);
        }

        // THIS IS JUST TO PRINT THE VECTOR OF DEVIATES FOR INDELS
        for (int x = 0; x < exp_deviates_vtr_ind.size(); x++) {
            std::cout << "i:" << exp_deviates_vtr_ind[x] << ' ';
        }

    }

    //printing out this stuff to check its working
    //std::cout << min_mut << "\n" << min_ind << "\n" << smallest_num << "\n" << position <<  "\n";
    
    std::cout << "\n" << "after mutateSeqEXP:\t" << simulated_protein << "\n" << "\n" ;
    return simulated_protein;
}

double getNormalDev(double mu, double stdev) {
    Normaldev mynorm(mu, stdev, myran.int64());
    double dev = mynorm.dev();
    //std::cout << dev << "\n";
    return dev;
}

//int main() {
//    std::string x = "MKNHCHKISAKHHHHAM";
//    mutateSeqExp(x); 
//}

////////// Changed here down - Brian ////////////
std::string mutateSeqExpBG(std::string simulated_protein, float mutation_rate, float indel_rate, float iterates){

    struct element {
        char aa;
        double mut;
        double ind;
        int length;
    };
    // Setting up the vectors
    vector<element> protein; 
    protein.push_back(element()); // create an element
    char aminoAcids[20] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E',
        'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };
                                      
    std::cout << "before mutateseqEXP:\t" << simulated_protein << "\n"; // Initially printing the non-mutated strin.
    int len = simulated_protein.length();
    for (int i = 0; i < simulated_protein.length(); i++) {
        protein.push_back(element()); // create an element
        protein[i].aa=simulated_protein[i];
    }

    // Mutation and indel rate set here now - dont need to set here
    //float mutation_rate = 0.14;
    //float indel_rate = 0.14;

    for (int i = 0; i < len; i++) {
        float beta1 = mutation_rate ; 
        Expondev myexp(beta1,myran.int64());
        double deviate = myexp.dev(); // choose exp_deviate(mean of beta)
        protein[i].mut=deviate;
    }
    // Traverse the string
    for (int i = 0; i < len; i++) {
        int counter = 1 ;
        //Code to scan back and forth to find repeats
        if (protein[i].aa != protein[i+1].aa && protein[i].aa != protein[i-1].aa) {
            float beta2 = indel_rate ; // the length if no repeats is 1
            Expondev myexp(beta2,myran.int64());
            double deviate = myexp.dev();
            protein[i].ind=deviate ; 
            protein[i].length=1; 
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
            float beta3 = indel_rate * counter ;
            Expondev myexp(beta3,myran.int64());
            double deviate = myexp.dev();
            protein[i].ind=deviate;
        }
    }

    for(int iter=0; iter<0.5*len; iter++) { // throw down 0.5 mut/site
        //selecting the lowest deviate from both vectors
        double minExpDev=protein[0].mut;
        int    minPosition=0;
        int    minType=0; // type=0 for mut, 1 for ins, 2 for del
        for (int i = 0; i < len; i++) {
            if(minExpDev > protein[i].mut) {
                minExpDev=protein[i].mut;
                minPosition=i; minType=0;
            }
            if(minExpDev > protein[i].ind) {
                minExpDev=protein[i].ind;
                minPosition=i; minType=1;
                if(myran.int64() < 0.5) { minType=2; }
            }
        }
        if(minType==0) { // put in a base change
            protein[minPosition].aa=aminoAcids[myran.int64() % numAA];
            float beta1 = mutation_rate ; 
            Expondev myexp(beta1,myran.int64());
            double deviate = myexp.dev();
            protein[minPosition].mut=deviate;
        } else { // put in an indel
            if(minType==1) { // use insertion
                if(len < 1.5*len) {
                    len++;
                    protein.push_back(element()); // create an element
                    for(int k=minPosition+1; k<len; k++) {
                        protein[k]=protein[k-1];
                    }
                    float beta1 = mutation_rate ; 
                    Expondev myexp(beta1,myran.int64());
                    double deviate = myexp.dev();
                    protein[minPosition].mut=deviate;
                }
            } else { // use deletion
                if(len > 0.5*len) {
                    len--;
                    for(int k=minPosition; k<len-1; k++) {
                        protein[k]=protein[k+1];
                    }
                }
            }
        }
        // check what is around the change
        // need to look down one, the site and up one
        for(int l=-1; l<2; l++) { // l=-1,0,1
            int i=minPosition+l;
            if (protein[i].aa == protein[i+1].aa || protein[i].aa == protein[i-1].aa) { // part of a repeat
                int x = i-1;
                //Look backward
                while (protein[i].aa == protein[x].aa && x >= 0) { x--; } 
                int start=x;
                //count forward
                x=start; int counter=1;
                while (protein[start].aa == protein[x].aa && x < len) { 
                    counter++; x++; 
                }
                int end=x;
                float beta3 = indel_rate * counter ;
                for(int j=start; j<=end; j++) {
                    Expondev myexp(beta3,myran.int64());
                    double deviate = myexp.dev();
                    protein[j].ind=deviate;
                    protein[j].length=counter;
                }
            } else { // not part of a repeat
                float beta3 = indel_rate;
                Expondev myexp(beta3,myran.int64());
                double deviate = myexp.dev();
                protein[i].ind=deviate;
                protein[i].length=1;
            }
        }
    }

    for(int i=0; i<len; i++) { simulated_protein[i]=protein[i].aa; }
    std::cout << "\n" << "after mutateSeqEXP:\t" << simulated_protein << "\n" << "\n" ;
    return simulated_protein;
}

