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

L=48

A=1.0
D=0.125
B=0.0
scale=1

FIT_LIMIT = 1e-4

term2s(x) = A*(2**(-D))**(D*exp(-(D+1)*x*scale) + D*D*(D+1)*exp(-x*scale*(D+3)))
term2c(x) = B*(2**(D-2))*((2-D)*exp(-((2-D)+1)*x*scale) + (2-D)*(2-D)*((2-D)+1)*exp(-x*scale*((2-D)+3)))


fit [10:48] term2s(x) + term2c(x) 'correlators_FTl1.dat' using 1:2:3 via A,D,B,scale

plot [0:48] "correlators_FTl1.dat" using 1:2:3 with yerrorbars  t "", term2s(x) + term2c(x) t ""
