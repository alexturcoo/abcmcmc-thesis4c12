#include <iostream>
#include <string.h>
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
        protein += aminoAcids[rand() % numAA];} //this rand() % 20 means in the range 0-19
    
    //std::cout << protein << "\n" << "\n" ;
    return protein;
}

//int main() {
//    createSeq(10);
//}
