set title "ABC-MCMC"
set xlabel "index"
set ylabel "Parameter Value"
set grid

plot "parameters.txt" u 1:2 w l, "parameters.txt" u 1:3 w l, "parameters.txt" u 1:3 w l
