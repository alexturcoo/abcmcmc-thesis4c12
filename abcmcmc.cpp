#include "mutations_2_BG_vecs.cpp"
#include "observed_protein.cpp"
#include "simulated_protein.cpp"
#include "distance2.cpp"
#include "vecAvg.cpp"
#include "ttest.cpp"
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

    double mut_rate_arr[10000];
    double ind_rate_arr[10000];
    double index[10000];

    //First for loop is for the number of simulations
    for (int i = 0; i < 1000; i++) {
    
    //setting this in the loop
    std::vector<vector<double>> vec_of_vecs; //for the current state
    std::vector<vector<double>> vec_of_vecs2; //for proposed state
    std::vector<double> distances;
    //importing the vector of summary statistic for observed protein SRP40 (Saccharomyces)
    //AND NORMALIZING IT, IT DOES NOT CHANGE
    std::vector<double> obs_prot_vtr = og_protein();
    obs_prot_vtr = normalize_vector(obs_prot_vtr);

        // 0. Proposing new parameter values
        double new_mut_rate = abs(getNormalDev2(mutation_rate, 1.0)); //+ mutation_rate;
        double new_indel_rate = abs(getNormalDev2(indel_rate, 1.0)); //+ indel_rate;
        std::cout << "Mutation rate: " << new_mut_rate << "\n";
        std::cout << "Indel Rate: " << new_indel_rate << "\n";

        // 1. We need to generate the random protein sequence
        // For loop here is to generate 10 vectors of summary
        // statistics for each parameter and get the average of all
        for (int k = 0; k<10; k++){
        std::string simulated_protein = createSeq(400); //For the current state
        std::string simulated_protein2 = createSeq(400); //For the proposed state
        
            // 2. Next we need to mutate the simulated protein
            // Going to try and mutate over 2 gens
            // this for loop is just to mutate the protein
            // Only going thru 1 mutation process right now - Feb 17
            // Might need to stick with a lower number of mutations,
            // around 200, 1000 giving me errors and idk why...
            for (int j = 0; j < 2; j++){
                std::string mutated_protein = mutateSeqExpBG(simulated_protein, mutation_rate, indel_rate);
                simulated_protein = mutated_protein;
            }
            // Also mutating a protein under the newly proposed
            // parameter values with this loop
            for (int k = 0; k < 2; k++) {
                std::string mutated_protein2 = mutateSeqExpBG(simulated_protein2, new_mut_rate, new_indel_rate);
                simulated_protein2 = mutated_protein2;
            }
        
            // 3. Next we will get the average of 10 vectors of summary statistics
            // For both original parameters and newly proposed
            // parameters
            std::vector<double> sim_prot_vtr = sim_protein(simulated_protein);
            vec_of_vecs.push_back(sim_prot_vtr);

            std::vector<double> sim_prot_vtr_2 = sim_protein(simulated_protein2);
            vec_of_vecs2.push_back(sim_prot_vtr_2);

            // Getting the 10 distances in case we need t-test
            std::vector<double> vec_normal = normalize_vector(sim_prot_vtr_2);
            double dist = vectors_distance2(vec_normal, obs_prot_vtr);
            distances.push_back(dist);
        }

        // TO TEST THE OUTPUT/PRINT THE VECTOR
        //for (int x = 0; x < vec_of_vecs.size(); x++) {
        //    for (int y = 0; y < vec_of_vecs[i].size(); y++){
        //        std::cout << vec_of_vecs[x][y] << " " << "\n";
        //    }
        //}

        // Getting the average vector of the 10 simulated vectors
        // for both original paramters and newly proposed parameters
        std::vector<double> sim_prot_vtravg = vectors_average(vec_of_vecs);
        std::vector<double> sim_prot_vtravg2 = vectors_average(vec_of_vecs2);

        // ADDED STEP - LETS NORMALIZE BOTH THE SIMULATED AND
        // OBSERVED VECTORS
        sim_prot_vtravg = normalize_vector(sim_prot_vtravg);
        sim_prot_vtravg2 = normalize_vector(sim_prot_vtravg2);

        // 5. Lets try finding the distance between the vectors
        double distance_current = vectors_distance2(sim_prot_vtravg, obs_prot_vtr);
        double distance_new = vectors_distance2(sim_prot_vtravg2, obs_prot_vtr);

        // 6. Does the distance satisfy the conditions?
        // If it does, add new parameter to corresponding list
        // otherwise stay at same value
        // Updated the code now to calculate the current state of
        // the model and the proposed state of the model to compare
        // those distances.
        if (distance_new < distance_current) {
            mutation_rate = new_mut_rate;
            indel_rate = new_indel_rate;
            mut_rate_arr[i] = mutation_rate;
            ind_rate_arr[i] = indel_rate;
            index[i] = i;
            std::cout << "ACCPETED" << "\n" << "\n";
        } else {//Do a ttest
            double p_value = tTest(distances, distance_current);
            double random_number = myran.doub(); //Random number 0-1
            if (random_number < p_value) {
                std::cout << "ACCEPTED pval" << "\n" << "\n";
                mutation_rate = new_mut_rate;
                indel_rate = new_indel_rate;
                mut_rate_arr[i] = mutation_rate;
                ind_rate_arr[i] = indel_rate;
                index[i] = i;
            } else {
                std::cout << "NOT ACCEPTED pval" << "\n" << "\n";
                mut_rate_arr[i] = mutation_rate;
                ind_rate_arr[i] = indel_rate;
                index[i] = i;
                continue;
            }
        }
  }

    std::ofstream myfile("parameters.txt"); //Create and open txt file
    // Printing the arrays of parameter values
    for (int b = 0; b<1000; b++) {
        myfile << index[b] << '\t' << mut_rate_arr[b] << '\t' << ind_rate_arr[b] << '\n';
    }
}

