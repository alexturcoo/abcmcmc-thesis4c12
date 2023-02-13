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


/////// FUNCTION TO GENERATE NORMAL DEVIATES //////
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

////////// Changed here down - Brian - Feb 12 Xander made Changes, trying to understand the new stuff ////////////
std::string mutateSeqExpBG(std::string simulated_protein, float mutation_rate, float indel_rate){ //took out iterates as a parameter

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
    // Traverse the string and generate deviates
    for (int i = 0; i < len; i++) {
        int counter = 1 ;
        //Code to scan back and forth to find repeats - this i-1 will give -1 on the first iteration though - check this in while loop too for scanning backwards
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
            float beta3 = indel_rate * counter ; // trying to see if the # of repeats plays a role, multiply by counter
            Expondev myexp(beta3,myran.int64());
            double deviate = myexp.dev();
            protein[i].ind=deviate;
        }
    }

    // I BELIEVE THIS PART IS CAUSING ERRORS RIGHT NOW - FEB 13
    // Working with the vector of structs, why we doing 0.5*len
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
        // So here were checking if the point mutation, insertion,
        // or deletion affected the landscape of the DNA sequence.
        // Did it create a repeat, inturrupt a repeat, etc.
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

///// TRYING TO DEBUG THIS FILE WITH THIS MAIN FUNCTION /////
int main() {
    for (int j = 0; j < 5; j++){
        std::string simulated_protein = createSeq(200);
        double mutation_rate = 0.14;
        double indel_rate = 0.14;
        std::string mutated_protein = mutateSeqExpBG(simulated_protein, mutation_rate, indel_rate);
        simulated_protein = mutated_protein;
    }        
}
