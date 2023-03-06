set title "ABC-MCMC"
set xlabel "index"
set ylabel "Parameter Value"
set grid

plot "accepted_parameters.txt" u 1:2 w l, "accepted_parameters.txt" u 1:3 w l
