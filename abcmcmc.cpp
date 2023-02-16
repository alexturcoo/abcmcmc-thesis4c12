#include "mutations_2_BG_vecs.cpp"
#include "observed_protein.cpp"
#include "simulated_protein.cpp"
#include "distance2.cpp"
#include "vecAvg.cpp"
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
    double mutation_rate = 0.14;
    double indel_rate = 0.14;
    double threshold = 0.005;

    double mut_rate_arr[1000];
    double ind_rate_arr[1000];
    double index[1000];

    //First for loop is for the number of simulations
    for (int i = 0; i < 1; i++) {
    
    //setting this in the loop
    std::vector<vector<double>> vec_of_vecs;

        // 0. Proposing new parameter values
        double new_mut_rate = getNormalDev(0.0, 1.0) + mutation_rate;
        double new_indel_rate = getNormalDev(0.0, 1.0) + indel_rate;
        std::cout << "Mutation rate: " << new_mut_rate << "\n";
        std::cout << "Indel Rate: " << new_indel_rate << "\n";

        // 1. We need to generate the random protein sequence
        // For loop here is to generate 10 vectors of summary
        // statistics for each parameter and get the average of all
        for (int k = 0; k<10; k++){
        std::string simulated_protein = createSeq(10);
        
            // 2. Next we need to mutate the simulated protein
            // Going to try and mutate over 2 gens
            // this for loop is just to mutate the protein
            for (int j = 0; j < 2; j++){
                std::string mutated_protein = mutateSeqExpBG(simulated_protein, mutation_rate, indel_rate);
                simulated_protein = mutated_protein;
            }
        
            // 3. Next we will get the average of 10 vectors of summary statistics
            //std::vector<double> sim_prot_vtr = sim_protein(simulated_protein);
            //vec_of_vecs.push_back(sim_prot_vtr);
        }

        // TO TEST THE OUTPUT/PRINT THE VECTOR
        //for (int x = 0; x < vec_of_vecs.size(); x++) {
        //    for (int y = 0; y < vec_of_vecs[i].size(); y++){
        //        std::cout << vec_of_vecs[x][y] << " " << "\n";
        //    }
        //}

        /*std::vector<double> sim_prot_vtravg = vectors_average(vec_of_vecs);

        // 4. Lets try importing the vector of summary statistic for observed protein SRP40 (Saccharomyces)
        std::vector<double> obs_prot_vtr = og_protein();

        // ADDED STEP - LETS NORMALIZE BOTH THE SIMULATED AND
        // OBSERVED VECTORS
        sim_prot_vtravg = normalize_vector(sim_prot_vtravg);
        obs_prot_vtr = normalize_vector(obs_prot_vtr);

        // 5. Lets try finding the distance between the vectors
        double distance = vectors_distance2(sim_prot_vtravg, obs_prot_vtr);

        // 6. Does the distance satisfy the conditions?
        // If it does, add new parameter to corresponding list
        // otherwise stay at same value
        if (distance < threshold) {
            mutation_rate = new_mut_rate;
            indel_rate = new_indel_rate;
            mut_rate_arr[i] = mutation_rate;
            ind_rate_arr[i] = indel_rate;
            index[i] = i;
            std::cout << "ACCPETED" << "\n" << "\n";
        } else {
            std::cout << "NOT ACCEPTED" << "\n" << "\n";
            mut_rate_arr[i] = mutation_rate;
            ind_rate_arr[i] = indel_rate;
            index[i] = i;
            continue;
        }*/

  }

    /*std::ofstream myfile("parameters.txt"); //Create and open txt file
    // Printing the arrays of parameter values
    for (int b = 0; b<10; b++) {
        myfile << index[b] << '\t' << mut_rate_arr[b] << '\t' << ind_rate_arr[b] << '\n';
    }*/
}

