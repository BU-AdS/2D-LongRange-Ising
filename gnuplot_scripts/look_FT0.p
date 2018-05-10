set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
set yrange       [0.0:0.3]
set terminal pdf
set output "FTl0_plot.pdf"
#set terminal dumb 120 40

set xlabel "t" font "courier,t"
set ylabel "Correlation function" font "courier,24"

L=48

A=0.5
D=0.125
B=1.0
scale=2*pi/(L/2)

FIT_LIMIT = 1e-12

term1s(x) = A*(exp(-D*L*scale)*cosh(D*(x-L)*scale) + D*D*exp(-(D+2)*L*scale)*cosh((D+2)*(x-L)*scale))
term1c(x) = B*(exp(-(2-D)*L*scale)*cosh((2-D)*(x-L)*scale) + (2-D)*(2-D)*exp(-((2-D)+2)*L*scale)*cosh(((2-D)+2)*(x-L)*scale))

fit [2:L] term1s(x) + term1c(x) 'correlators_FTl0.dat' using 1:2:3 via A,D,B

set label sprintf("C_{/Symbol c}/C_{/Symbol s} = %f\n\n{/Symbol D} = %f\n\na = %f", B/A, D, scale) at 25,0.25

plot [0:L] "correlators_FTl0.dat" using 1:2:3 with yerrorbars  t "", term1s(x) + term1c(x) t ""
