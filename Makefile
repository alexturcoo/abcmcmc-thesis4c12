abcmcmc: abcmcmc.o
	c++ -o abcmcmc abcmcmc.o -lm

abcmcmc.o: abcmcmc.cpp
	c++ -c abcmcmc.cpp

createseq: createseq.o
	c++ -o createseq createseq.o -lm

createseq.o: createseq.cpp
	c++ -c createseq.cpp

mutations: mutations.o
	c++ -o mutations mutations.o -lm

mutations.o: mutations.cpp
	c++ -c mutations.cpp

observed_protein: observed_protein.o
	c++ -o observed_protein observed_protein.o -lm

observed_protein.o: observed_protein.cpp
	c++ -c observed_protein.cpp
