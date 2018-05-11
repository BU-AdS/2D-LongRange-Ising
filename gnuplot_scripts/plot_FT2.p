set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
set yrange       [0.0:0.06]
set terminal pdf
set output "FTl2plot.pdf"
#set terminal dumb 120 40

set xlabel "t" font "courier,t"
set ylabel "Correlation function" font "courier,24"

L=48

A=1.0
D=0.125
B=1.0
scale=2*pi/(L/2)

FIT_LIMIT = 1e-6

term3s(x) = A*(D*(D+1)*exp(-(D+2)*L*scale)*cosh((D+2)*(x-L)*scale) + (1.0/3.0)*D*D*(D+1)*(D+2)*exp(-(D+4)*L*scale)*cosh((D+4)*(x-L)*scale))
term3c(x) = B*((2-D)*((2-D)+1)*exp(-((2-D)+2)*L*scale)*cosh(scale*((2-D)+2)*(x-L)) + (1.0/3.0)*(2-D)*(2-D)*((2-D)+1)*((2-D)+2)*exp(-((2-D)+4)*L*scale)*cosh(((2-D)+4)*(x-L)*scale))

fit [3:L] term3s(x) + term3c(x) 'correlators_FTl2.dat' using 1:2:3 via A,B,D,scale

plot [0:L] "correlators_FTl2.dat" using 1:2:3 with yerrorbars  t "", term3s(x) + term3c(x) t ""
