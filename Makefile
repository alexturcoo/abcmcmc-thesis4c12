abcmcmc: abcmcmc.o
	c++ -o abcmcmc abcmcmc.o -lm

abcmcmc.o: abcmcmc.cpp
	c++ -c abcmcmc.cpp

createseq: createseq.o
	c++ -o createseq createseq.o -lm

createseq.o: createseq.cpp
	c++ -c createseq.cpp

mutations_2: mutations_2.o
	c++ -o mutations_2 mutations_2.o -lm

mutations_2.o: mutations_2.cpp
	c++ -c mutations_2.cpp

observed_protein: observed_protein.o
	c++ -o observed_protein observed_protein.o -lm

observed_protein.o: observed_protein.cpp
	c++ -c observed_protein.cpp

distance: distance.o
	c++ -o distance distance.o -lm

distance.o: distance.cpp
	c++ -c distance.cpp

distance2: distance2.o
	c++ -o distance2 distance2.o -lm

distance2.o: distance2.cpp
	c++ -c distance2.cpp

simulated_protein: simulated_protein.o
	c++ -o simulated_protein simulated_protein.o -lm

simulated_protein.o: simulated_protein.cpp
	c++ -c simulated_protein.cpp

getindex: getindex.o
	c++ -o getindex getindex.o -lm

getindex.o: getindex.cpp
	c++ -c getindex.cpp
