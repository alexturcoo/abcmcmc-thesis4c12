#include <iostream>
#include <fstream>
#include <string>


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

int main() {
    og_protein();
}
