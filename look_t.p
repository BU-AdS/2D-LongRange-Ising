set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
#set yrange       [-2:0]
#set terminal postscript color enhanced
set terminal dumb 120 40

coeff=0.4343

M=1
C=1

lin(x) = -M*x + C

fit [0:1.5] lin(x) 'correlators_t.dat' using (log($1)):(log($2)):(coeff*(($3)/($2))) via M,C

plot "correlators_t.dat" using (log($1)):(log($2)):(coeff*(($3)/($2))) with yerrorbars, lin(x)
