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

    // Setting initial parameters
    double mutation_rate = 10.00;
    double indel_rate = 10.00;
    double threshold = 0.5;

    // For the new proposals
    double mean = 0.00;
    double stdev = 1.00;

    std::vector<double> mut_rate_vtr;
    std::vector<double> ind_rate_vtr;

    //First for loop is for the number of simulations
    for (int i = 0; i < 1000000; i++) {
        // 0. Proposing new parameter values
        double new_mut_rate = getNormalDev(mean, stdev) + mutation_rate;
        double new_indel_rate = getNormalDev(mean, stdev) + indel_rate;

        // 1. We need to generate the random protein sequence
        std::string simulated_protein = createSeq(400);

        // 2. Next we need to mutate the simulated protein
        // Going to try and mutate over 2 gens
        // this for loop is just to mutate the protein
        for (int j = 0; j < 4000; j++){
            std::string mutated_protein = mutateSeqExp(simulated_protein);
            simulated_protein = mutated_protein;
        }

        // 3. Next we will get the vector of summary statistic for simulated protein
        std::vector<double> sim_prot_vtr = sim_protein(simulated_protein);

        // 4. Lets try importing the vector of summary statistic for observed protein
        std::vector<double> obs_prot_vtr = og_protein();

        // 5. Lets try finding the distance between the vectors
        double distance = vectors_distance2(sim_prot_vtr, obs_prot_vtr);

        // 6. Does the distance satisfy the conditions?
        if (distance < threshold) {
            mutation_rate = new_mut_rate;
            indel_rate = new_indel_rate;
            mut_rate_vtr.push_back(mutation_rate);
            ind_rate_vtr.push_back(indel_rate);
        } else {
            continue;
        }
    }
}

