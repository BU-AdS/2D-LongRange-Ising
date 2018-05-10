set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
set yrange       [0.0:0.2]
set terminal pdf
set output "FTl1_plot.pdf"
#set terminal dumb 120 40

set xlabel "t" font "courier,t"
set ylabel "Correlation function" font "courier,24"

L=48

A=1.0
D=0.125
B=1.0
scale=2*pi/(L/2)

FIT_LIMIT = 1e-4

term2s(x) = A*D*exp(-(D+1)*L*scale)*cosh((D+1)*(x-L)*scale) + D*D*(D+1)*exp(-(D+3)*L*scale)*cosh((D+3)*(x-L)*scale) )
term2c(x) = B*((2-D)*exp(-((2-D)+1)*L*scale)*cosh(scale*((2-D)+1)*(x-L)) + (2-D)*(2-D)*((2-D)+1)*exp(-((2-D)+3)*L*scale)*cosh(((2-D)+3)*(x-L)*scale))

fit [1:L] term2s(x) + term2c(x) 'correlators_FTl1.dat' using 1:2:3 via A,D,B,scale

plot [0:L] "correlators_FTl1.dat" using 1:2:3 with yerrorbars  t "", term2s(x) + term2c(x) t ""
