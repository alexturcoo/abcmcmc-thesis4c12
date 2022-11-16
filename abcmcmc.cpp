#include "createseq.cpp"
#include "mutations_2.cpp"
#include "observed_protein.cpp"
#include "simulated_protein.cpp"
#include "distance2.cpp"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <bits/stdc++.h>
#include <fstream> //The fstream library allows us to work with files
#include <vector> // Working with vectors is nicer than arrays, more c++

//////////////////////////////////////////////
// THE MAIN FUNCTION - WHERE EXECUTION BEGINS/
//////////////////////////////////////////////

int main() {

    // 1. We need to generate the random protein sequence
    std::string simulated_protein = createSeq(15);
    // 2. Next we need to mutate the simulated protein
    // Going to try and mutate over 2 gens
    for (int i = 0; i < 10; i++){
        std::string mutated_protein = mutateSeqExp(simulated_protein);
        simulated_protein = mutated_protein;
    }
    // 3. Next we will get the vector of summary statistic for simulated protein
    std::vector<double> sim_prot_vtr = sim_protein(simulated_protein);
    // 4. Lets try importing the vector of summary statistic for observed protein
    std::vector<double> obs_prot_vtr = og_protein();
    //5. Lets try finding the distance between the vectors
    double distance = vectors_distance2(sim_prot_vtr, obs_prot_vtr);
}

