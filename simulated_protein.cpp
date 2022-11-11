#include <iostream>
#include <fstream>
#include <string>
#include <vector>

////////////////////////////////////////////////////////////////////
// SECOND FUNCTION READS IN SRP40 PROTEIN AND CREATES A VECTOR WITH/
// CHARACTERISTICS - THIS IS EITHER S or D - THE ORIGINAL DATA     /
////////////////////////////////////////////////////////////////////

std::vector<double> sim_protein(std::string sim_prot){

    std::vector<double> sim_ss_vec; //Making the simulated summary statistic vector

    // Finding the LENGTH of protein sequence
    int length = sim_prot.length();
    sim_ss_vec.push_back(length);

    // Creating fasta file of simulated protein string to use in segA
    std::ofstream MyFile("simulated_protein.fasta"); //Create and open a fasta file
    MyFile << ">simulation\n" << sim_prot; //Write to the file
    MyFile.close(); //Close the file

    // Finding the number of LCRs
    system("segA simulated_protein.fasta 15 1.9 2.2 -l | grep '>' | wc -l >> numlcr_sim.txt"); //Write the results of counting the LCRs to a txt file
    std::ifstream myreadfile3("numlcr_sim.txt"); //Declaring fstream variable and connecting it to a stream object by opening the file
    double numlcr_sim; //Declare int variable which holds the contents of the 'numlcr.txt' file
    myreadfile3 >> numlcr_sim; //Pipe file's content into stream
    sim_ss_vec.push_back(numlcr_sim); // using push_back to append to vector
    
    //This if statement is to remove the division by zero error in
    //the awk command below
    if (numlcr_sim == 0) {
        sim_ss_vec.push_back(0);
        //Removing files so they are empty after
        system("rm numlcr_sim.txt");
        system("rm simulated_protein.fasta");
    } else {

    // Finding the average entropy of the LCRs dont know why I have
    // to use egrep for the extraction of numbers part in the system
    // call
    system("segA simulated_protein.fasta 15 1.9 2.2 -l | grep 'complexity' | awk '{print $2}' | egrep -o '([0.0-9.0]+)' | awk '{sum += $1} END {print mean=sum/NR}' >> avg_entropy_lcrs_sim.txt");
    std::ifstream myreadfile4("avg_entropy_lcrs_sim.txt");
    double avgentropylcrs_sim; 
    myreadfile4 >> avgentropylcrs_sim;
    sim_ss_vec.push_back(avgentropylcrs_sim);
    //Removing files so they are empty after
    system("rm avg_entropy_lcrs_sim.txt");
    system("rm numlcr_sim.txt");
    system("rm simulated_protein.fasta");
    }

    // TO TEST THE OUTPUT/PRINT THE VECTOR
    for (int x = 0; x < sim_ss_vec.size(); x++) {
        std::cout << sim_ss_vec[x] << ' ' << "\n";
    }

    return sim_ss_vec;
}

//int main() {

    //sim_protein("MASKKIKVDEVPKLSVKEKEIEEKSSSSSSSSSSSSSSSSSSSSSSSSSGESSSSSSSSSSSSSSDSSDSSDSESSSSSSSSSSSSSSSSDSESSSESDSSSSGSSSSSSSSSDESSSESESEDETKKRARESDNEDAKETKKAKTEPESSSSSESSSSGSSSSSESESGSESDSDSSSSSSSSSDSESDSESDSQSSSSSSSSDSSSDSDSSSSDSSSDSDSSSSSSSSSSDSDSDSDSSSDSDSSGSSDSSSSSDSSSDESTSSDSSDSDSDSDSGSSSELETKEATADESKAEETPASSNESTPSASSSSSANKLNIPAGTDEIKEGQRKHFSRVDRSKINFEAWELTDNTYKGAAGTWGEKANEKLGRVRGKDFTKNKNKMKRGSYRGGSITLESGSYKFQD");
//    sim_protein(createSeq(400));
//    return 0;
//}
