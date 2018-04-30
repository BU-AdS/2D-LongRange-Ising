set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
set yrange       [0.0:0.06]
set terminal pdf
set output "plot.pdf"
#set terminal dumb 120 40

set xlabel "t" font "courier,t"
set ylabel "Correlation function" font "courier,24"

L=48

A=1.0
D=0.125
B=1.0
scale=1

FIT_LIMIT = 1e-4

term3s(x) = A*(2**(-D))**(D*(D+1)*exp(-x*(D+2)*scale) + (1.0/3.0)*D*D*(D+1)*(D+2)*exp(-x*(D+4)*scale))
term3c(x) = B*(2**(-(2-D)))**((2-D)*((2-D)+1)*exp(-x*((2-D)+2)*scale) + (1.0/3.0)*(2-D)*(2-D)*((2-D)+1)*((2-D)+2)*exp(-x*((2-D)+4)*scale))

fit [10:48] term3s(x) + term3c(x) 'correlators_FTl2.dat' using 1:2:3 via A,D,B,scale

plot [0:48] "correlators_FTl2.dat" using 1:2:3 with yerrorbars  t "", term3s(x) + term3c(x) t ""