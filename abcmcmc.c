#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include "myRandomFunctions.c"

// THIS IS A TEST TO SEE IF I CAN PUSH AND PULL FROM GITHUB
// PUT THINGS YOURE GONNA DEFINE HERE
// 21 amino acids because of 0 at the end which is not an amino acid
#define numAA 20

// FIRST FUNCTION TAKES AN INTEGER VALUE AND GENERATES A RANDOM
// AMINO ACID SEQUENCE OF THAT LENGTH

std::string createSeq(int n){

    char aminoAcids[numAA] = { 'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 
                                   'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T' };

    std::string protein = "";
    for (int i = 0; i < n; i++){
        protein = protein + aminoAcids[rand() % numAA];}
    
    std::cout << protein << "\n" << "\n" ;
   // printf("%s", protein.c_str());
    return protein;
}

// SECOND FUNCTION READS IN SRP40 PROTEIN AND CREATES A VECTOR WITH
// CHARACTERISTICS - THIS IS EITHER S or D - THE ORIGINAL DATA

int og_protein(){

    std::string srp40 = "MASKKIKVDEVPKLSVKEKEIEEKSSSSSSSSSSSSSSSSSSSSSSSSSGESSSSSSSSSSSSSSDSSDS\
SDSESSSSSSSSSSSSSSSSDSESSSESDSSSSGSSSSSSSSSDESSSESESEDETKKRARESDNEDAKE\
TKKAKTEPESSSSSESSSSGSSSSSESESGSESDSDSSSSSSSSSDSESDSESDSQSSSSSSSSDSSSDS\
DSSSSDSSSDSDSSSSSSSSSSDSDSDSDSSSDSDSSGSSDSSSSSDSSSDESTSSDSSDSDSDSDSGSS\
SELETKEATADESKAEETPASSNESTPSASSSSSANKLNIPAGTDEIKEGQRKHFSRVDRSKINFEAWEL\
TDNTYKGAAGTWGEKANEKLGRVRGKDFTKNKNKMKRGSYRGGSITLESGSYKFQD";

    int length = srp40.length();
    int num_lcr = 
    int avg_entropy_lcr =

    std::cout << srp40 << "\n" << "\n" << length << "\n";
    return 0;

}

// THE MAIN FUNCTION - WHERE EXECUTION BEGINS
int main() {
    
    srand(18);
    createSeq(100); // Creating the Protein
    og_protein();   // Generating data for SRP40 protein
    return 0;
}

