set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
set yrange       [0.0:0.2]
set terminal pdf
set output "plot.pdf"
#set terminal dumb 120 40

set xlabel "t" font "courier,t"
set ylabel "Correlation function" font "courier,24"

A=1.0
D=0.125
B=1.0
scale=2*pi/24

FIT_LIMIT = 1e-12

term1s(x) = A*(2**(-D))*(exp(-x*D*scale) + D*D*exp(-x*(D+2)*scale) + 0.25*D*D*(D+1)*(D+1)*exp(-x*(D+4)*scale))
term1c(x) = B*(2**(D-2))*(exp(-x*scale*(2-D)) + (2-D)*(2-D)*exp(-x*scale*((2-D)+2)) + 0.25*(2-D)*(2-D)*((2-D)+1)*((2-D)+1)*exp(-x*((2-D)+4)*scale))

fit [8:40] term1s(x) + term1c(x) 'correlators_FTl0.dat' using 1:2:3 via A,D,B

plot [0:48] "correlators_FTl0.dat" using 1:2:3 with yerrorbars  t "", term1s(x) + term1c(x) t ""
