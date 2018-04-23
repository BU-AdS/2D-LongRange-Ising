set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
set xrange       [0:10]
set yrange       [-0.1:0.5]
set terminal pdf
set output "SigCrit.pdf"
#set terminal dumb 120 40

A=1
B=1
C=-1.75

#fit(x) = A*sqrt(B*(x-C))
fit(x) = A*log(B*(x+C))

fit [2.5:12] fit(x) 'critical.dat' using 1:2 via B,A

level(x) = A*log(B*(5+C))

plot [0:12] "critical.dat" using 1:2, fit(x), level(x)
