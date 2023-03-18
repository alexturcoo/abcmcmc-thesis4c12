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
    double mutation_rate = 1.00;
    double indel_rate = 1.00;
    int num_simulations = 1;
    int num_mutations = 1;
    double mean_proposal = 0.0;
    double stddev_proposal = 1.0;

    double mut_rate_arr[10000]; //for accepted mutation rates
    double ind_rate_arr[10000]; //for accepted indel rates
    double index[10000]; // simulation iterations
    double distance_array[10000]; //for accepted distances
    
    double current_mut_rate_arr[10000]; //for current mutation rates
    double proposed_mut_rate_arr[10000]; //for proposed mutation rates
    double current_ind_rate_arr[10000]; //for current indel rates
    double proposed_ind_rate_arr[10000]; //for proposed indel rates
    double current_distances_arr[10000]; //for current distances avg of 10
    double proposed_distances_arr[10000]; //for proposed distances avg of 10
    double distances_ttest[10000]; //for each of the 10 distances for proposed parameters

    double acc_rej_rate[10000]; //acceptance-rejection rate for mutation
    
    double mean_ttest_arr[10000]; //storing means from ttest of 10 vectors
    double sttdev_ttest_arr[10000]; //storing stdev from 10 vectors
    double tstat_array[10000]; //storing test statistics
    double probability_ttest[10000]; //storing pvalues for ttest
    double accepted_rejected[10000]; //1 for accepted, 0 for rejected
    std::vector<vector<double>> vec_of_vecs3; //for getting vectors of 10 distances
    std::vector<vector<double>> vec_of_vecs4; //summary statistics for each of 10

    std::vector<vector<double>> vec_of_vecs_length;
    std::vector<vector<double>> vec_of_vecs_avgnum;
    std::vector<vector<double>> vec_of_vecs_entropy;

    double new_mut_rate = 0;
    double new_ind_rate = 0;
    double acceptance_counter = 0;

    //importing the vector of summary statistic for observed protein SRP40 (Saccharomyces)
    //AND NORMALIZING IT, IT DOES NOT CHANGE
    std::vector<double> obs_prot_vtr = og_protein();
    obs_prot_vtr = normalize_vector(obs_prot_vtr);

    //Generating current state protein sequence outside loop for
    //first time
    std::string simulated_protein = createSeq(400); //For the current state

    // GETTING THE FIRST CURRENT DISTANCE OUTSIDE THE MAIN LOOP
    // TRYING TO GET MUTATED SEQUENCE HERE SO WE CAN KEEP MUTATING
    // THE SAME SEQUENCE IN THE ALGORITHM
    for (int j = 0; j < num_mutations; j++){
        std::string mutated_protein = mutateSeqExpBG(simulated_protein, mutation_rate, indel_rate);
        simulated_protein = mutated_protein;
    }

    std::vector<double> sim_prot_vtr = sim_protein(simulated_protein);
    sim_prot_vtr = normalize_vector(sim_prot_vtr);
    double distance_current = vectors_distance2(sim_prot_vtr, obs_prot_vtr);
    std::cout << distance_current << '\n';

    //First for loop is for the number of simulations
    for (int i = 0; i < num_simulations; i++) {
    
        //setting this in the loop
        std::vector<vector<double>> vec_of_vecs; //for the current state
        std::vector<vector<double>> vec_of_vecs2; //for proposed state
        std::vector<double> distances; //holding vector of 10 distances per simulation iteration
        std::vector<double> length_vec; //holding length summary statistic
        std::vector<double> num_lcrs_vec; //holding number of LCRs summary statistic
        std::vector<double> avg_entropy_vec; //holding average entropy summary statistic

        current_mut_rate_arr[i] = mutation_rate; // adding current mut rate to array
        current_ind_rate_arr[i] = indel_rate; // adding current ind rate to array

        // 0. Proposing new parameter values - this if statement is
        // to hold one of the parameters constant on every other
        // iteration... in other words, holding one parameter
        // constant and changing the other on every other
        // iteration. - I changed this back to just generating new
        // ones each time.
        new_mut_rate = getNormalDev2(mean_proposal, stddev_proposal) + mutation_rate; //+ mutation_rate;
        std::cout << "Mutation rate: " << new_mut_rate << "\n";
        new_ind_rate = getNormalDev2(mean_proposal, stddev_proposal) + indel_rate; //+ indel_rate;i
        std::cout << "Indel Rate: " << new_ind_rate << "\n";

        if (i == 0) {
            new_mut_rate = mutation_rate;
            new_ind_rate = indel_rate;
        }

        proposed_mut_rate_arr[i] = new_mut_rate; // adding the proposed mut rate to array
        proposed_ind_rate_arr[i] = new_ind_rate; // adding the proposed ind rate to array

        std::string simulated_protein2 = createSeq(400); //for proposed state - going to use the same sequence and keep mutating it for all 10 iterates below

        // 1. We need to generate the random protein sequence
        // For loop here is to generate 10 vectors of summary
        // statistics for each parameter and get the average of all
        for (int k = 0; k<1000; k++){     
            // 2. Next we need to mutate the simulated protein
            // Going to try and mutate over 2 gens
            // this for loop is just to mutate the protein
            // Only going thru 1 mutation process right now - Feb 17
            // Might need to stick with a lower number of mutations,
            // around 200, 1000 giving me errors and idk why...
            // Also mutating a protein under the newly proposed
            // parameter values with this loop
            for (int h = 0; h < num_mutations; h++) {
                std::string mutated_protein2 = mutateSeqExpBG(simulated_protein2, new_mut_rate, new_ind_rate);
                simulated_protein2 = mutated_protein2;
            }
        
            // 3. Next we will get the average of 10 vectors of summary statistics
            // For newly proposed parameter values

            std::vector<double> sim_prot_vtr_2 = sim_protein(simulated_protein2);
            vec_of_vecs2.push_back(sim_prot_vtr_2);

            // Getting the 10 distances in case we need t-test
            std::vector<double> vec_normal = normalize_vector(sim_prot_vtr_2);
            double length = sim_prot_vtr_2[0];
            double numlcr = sim_prot_vtr_2[1];
            double avgent = sim_prot_vtr_2[2];
            length_vec.push_back(length);
            num_lcrs_vec.push_back(numlcr);
            avg_entropy_vec.push_back(avgent);
            double dist = vectors_distance2(vec_normal, obs_prot_vtr);
            distances.push_back(dist); //storing the 10 distances to get mean and stdev
            distances_ttest[k] = dist; //storing each of the 10 distances for output
            
        }

        vec_of_vecs3.push_back(distances); //pushing back all 10 distances
        vec_of_vecs_length.push_back(length_vec); //pushing back lengths
        vec_of_vecs_avgnum.push_back(num_lcrs_vec);
        vec_of_vecs_entropy.push_back(avg_entropy_vec);

        // TO TEST THE OUTPUT/PRINT THE VECTOR
        //for (int x = 0; x < vec_of_vecs.size(); x++) {
        //    for (int y = 0; y < vec_of_vecs[i].size(); y++){
        //        std::cout << vec_of_vecs[x][y] << " " << "\n";
        //    }
        //}

        // Getting the average vector of the 10 simulated vectors
        // for newly proposed parameters
        std::vector<double> sim_prot_vtravg2 = vectors_average(vec_of_vecs2);

        // ADDED STEP - LETS NORMALIZE THE VECTOR
        sim_prot_vtravg2 = normalize_vector(sim_prot_vtravg2);

        // 5. Lets try finding the distance between the vectors
        double distance_new = vectors_distance2(sim_prot_vtravg2, obs_prot_vtr);

        current_distances_arr[i] = distance_current; //storing current distances
        proposed_distances_arr[i] = distance_new; //storing proposed new distances

        // 6. Does the distance satisfy the conditions?
        // If it does, add new parameter to corresponding list
        // otherwise stay at same value
        // Updated the code now to calculate the current state of
        // the model and the proposed state of the model to compare
        // those distances.

        if (distance_new < distance_current && new_mut_rate > 0 && new_ind_rate > 0) {
            mutation_rate = new_mut_rate;
            indel_rate = new_ind_rate;
            mut_rate_arr[i] = mutation_rate;
            ind_rate_arr[i] = indel_rate;
            distance_array[i] = distance_new;
            index[i] = i;
            accepted_rejected[i] = 1;
            acceptance_counter = acceptance_counter + 1;
            acc_rej_rate[i] = acceptance_counter/i;
            distance_current = distance_new; //setting the new current distance to the new distance so we can compare a new distance to this accepted one
            std::cout << "ACCPETED" << "\n" << "\n";
        } else {//Do a ttest
            double standard_dev = stdDev(distances);
            double mean = Mean(distances);
            double tstat = abs((mean - distance_current) / (standard_dev / sqrt(distances.size())));
            double p_value = tTest(distances, distance_current);

            //Adding ttest values to output
            mean_ttest_arr[i] = mean; //storing means from ttest of 10 vectors
            sttdev_ttest_arr[i] = standard_dev; //storing stdev from 10 vectors
            tstat_array[i] = tstat; //storing test statistics
            probability_ttest[i] = p_value; //storing pvalues for ttest

            double random_number = myran.doub(); //Random number 0-1
            if (random_number < p_value && new_mut_rate > 0 && new_ind_rate > 0) {
                std::cout << "ACCEPTED pval" << "\n" << "\n";
                mutation_rate = new_mut_rate;
                indel_rate = new_ind_rate;
                mut_rate_arr[i] = mutation_rate;
                ind_rate_arr[i] = indel_rate;
                distance_array[i] = distance_new;
                index[i] = i;
                accepted_rejected[i] = 1;
                acceptance_counter = acceptance_counter + 1;
                acc_rej_rate[i] = acceptance_counter/i;
                distance_current = distance_new;
            } else {
                std::cout << p_value << "\n" << random_number << "\n";
                std::cout << "NOT ACCEPTED pval" << "\n" << "\n";
                std::cout << i << "\n";
                mut_rate_arr[i] = mutation_rate;
                ind_rate_arr[i] = indel_rate;
                distance_array[i] = distance_current;
                index[i] = i;
                accepted_rejected[i] = 0;
                acc_rej_rate[i] = acceptance_counter/i;
                continue;
            }
        }
  }

    std::ofstream myfile("accepted_parameters.txt"); //Create and open txt file
    // Printing the arrays of parameter values
    for (int b = 0; b<num_simulations; b++) {
        myfile << index[b] << '\t' << mut_rate_arr[b] << '\t' << ind_rate_arr[b] << '\t' << distance_array[b] << '\t' << '\t' <<  current_mut_rate_arr[b] << '\t' << current_ind_rate_arr[b] << '\t' << 
            proposed_mut_rate_arr[b] << '\t' << proposed_ind_rate_arr[b] << '\t' << '\t' << mean_ttest_arr[b] << '\t' <<  sttdev_ttest_arr[b] << '\t' <<  tstat_array[b] << '\t' << probability_ttest[b] << '\t'  
            << accepted_rejected[b] << '\t' << '\t' << acc_rej_rate[b] << '\n';
    }

    std::ofstream myfile2("10_distances.txt");
    for (int x = 0; x<vec_of_vecs3.size(); x++) {
        for (int y = 0; y<vec_of_vecs3[x].size(); y++) {
            myfile2 << index[x] << '\t' << '\t' << proposed_mut_rate_arr[x] << '\t' << proposed_ind_rate_arr[x] << '\t' << vec_of_vecs3[x][y] << '\t' << vec_of_vecs_length[x][y]
                << '\t' << vec_of_vecs_avgnum[x][y] << '\t' << vec_of_vecs_entropy[x][y] << '\n';
        }
    }
}
