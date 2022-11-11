#include "createseq.cpp"
#include "mutations.cpp"
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
    std::string mutated_protein = mutateSeqAA(simulated_protein);
    mutated_protein = mutateSeqExp(mutated_protein);
    // 3. Lets try getting the vector of ss for simulated proteini
    sim_protein(mutated_protein);
}

