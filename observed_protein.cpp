#include <iostream>
#include <fstream>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////
// SECOND FUNCTION READS IN SRP40 PROTEIN AND CREATES A VECTOR WITH/
// CHARACTERISTICS - THIS IS EITHER S or D - THE ORIGINAL DATA     /
////////////////////////////////////////////////////////////////////

std::vector<double> og_protein(){

    std::string srp40 = "MASKKIKVDEVPKLSVKEKEIEEKSSSSSSSSSSSSSSSSSSSSSSSSSGESSSSSSSSSSSSSSDSSDSSDSESSSSSSSSSSSSSSSSDSESSSESDSSSSGSSSSSSSSSDESSSESESEDETKKRARESDNEDAKETKKAKTEPESSSSSESSSSGSSSSSESESGSESDSDSSSSSSSSSDSESDSESDSQSSSSSSSSDSSSDSDSSSSDSSSDSDSSSSSSSSSSDSDSDSDSSSDSDSSGSSDSSSSSDSSSDESTSSDSSDSDSDSDSGSSSELETKEATADESKAEETPASSNESTPSASSSSSANKLNIPAGTDEIKEGQRKHFSRVDRSKINFEAWELTDNTYKGAAGTWGEKANEKLGRVRGKDFTKNKNKMKRGSYRGGSITLESGSYKFQD";


    // Creating fasta file of protein string to use in segA
    std::ofstream MyFile("srp40.fasta"); //Create and open a fasta file
    MyFile << ">srp40_saccharomyces_cerevisiae\n" << srp40; //Write to the file
    MyFile.close(); //Close the file
    
    std::vector<double> obs_ss_vec; //Making the observed summary statistic vector

    // Finding the LENGTH of protein sequence
    int length = srp40.length();
    obs_ss_vec.push_back(length);

    // Finding the number of LCRs
    system("segA srp40.fasta 15 1.9 2.2 -l | grep '>' | wc -l >> numlcr.txt"); //Write the results of counting the LCRs to a txt file
    std::ifstream myreadfile("numlcr.txt"); //Declaring fstream variable and connecting it to a stream object by opening the file
    int num_lcr; //Declare int variable which holds the contents of the 'numlcr.txt' file
    myreadfile >> num_lcr; //Pipe file's content into stream
    obs_ss_vec.push_back(num_lcr); // using push_back to append to vector
    
    //This if statement is to remove the division by zero error in
    //the awk command below
    if (num_lcr == 0) {
        obs_ss_vec.push_back(0);
        //Removing files so they are empty after
        system("rm numlcr.txt");
        system("rm srp40.fasta");
    } else {

    // Finding the average entropy of the LCRs dont know why I have
    // to use egrep for the extraction of numbers part in the system
    // call
    system("segA srp40.fasta 15 1.9 2.2 -l | grep 'complexity' | awk '{print $2}' | egrep -o '([0.0-9.0]+)' | awk '{sum += $1} END {print mean=sum/NR}' >> avg_entropy_lcrs.txt");
    std::ifstream myreadfile2("avg_entropy_lcrs.txt");
    float avgentropylcrs; 
    myreadfile2 >> avgentropylcrs;
    obs_ss_vec.push_back(avgentropylcrs);

    // Removing the files so thy are empty for the next loop 
    system("rm avg_entropy_lcrs.txt");
    system("rm numlcr.txt");
    system("rm srp40.fasta");
    }

    // TO TEST THE OUTPUT/PRINT THE VECTOR
    for (int x = 0; x < obs_ss_vec.size(); x++) {
        std::cout << obs_ss_vec[x] << ' ' << "\n";
    }

    return obs_ss_vec;
}

//int main() {
//    og_protein();
//    return 0;
//}
