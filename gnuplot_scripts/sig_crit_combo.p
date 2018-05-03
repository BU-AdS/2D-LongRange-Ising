set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
set xrange       [0:10]
set yrange       [0.0:0.5]
set terminal pdf
set output "SigCrit.pdf"
#set terminal dumb 120 40

set xlabel "sigma" font "courier,24"
set ylabel "T crit" font "courier,24"

A=1
B=1
C=1
E=1
F=1
G=1


m=1
c=0.02
n=1
d=0.02

fit1(x) = m*x+c
fit2(x) = A*log(B*(x-C))

fit3(x) = n*x+d
fit4(x) = E*log(F*(x-G))
		
fit [1.1:1.75] fit1(x) 'critical_16.dat' using 1:2 via m,c
fit [1.75:10] fit2(x) 'critical_16.dat' using 1:2 via A,B,C

fit [1.1:1.75] fit3(x) 'critical_32.dat' using 1:2 via n,d
fit [1.75:10] fit4(x) 'critical_32.dat' using 1:2 via E,F,G

plot [0:11] "critical_16.dat" using 1:2 t "Critical line (16)", "critical_32.dat" using 1:2 t "Critical line (32)", fit1(x) t "16 Linear (sigma < 1.75)", fit2(x) t "16 logarithmic (sigma > 1.75)", fit3(x) t "32 Linear (sigma < 1.75)", fit4(x) t "32 logarithmic (sigma > 1.75)"
