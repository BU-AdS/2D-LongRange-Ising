set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
set yrange       [-0.1:1.0]
set terminal pdf
set output "plot.pdf"
#set terminal dumb 120 40

set xlabel "t" font "courier,t"
set ylabel "Correlation function" font "courier,24"

L=48

A=1.9
B=1.0
COS_T=cos(0*2*pi/(L/2))
scale=2*pi/(L/2)
a=1
b=2
c=3
f=4
g=5
h=6

D=0.125
d=1.875

FIT_LIMIT = 1e-12

term1s(x) = A*(2**(-D))*exp(-D*L*scale)*cosh(D*(x-L)*scale) 
term1c(x) = B*(2**(-(2-D)))*exp(-(2-D)*L*scale)*cosh((2-D)*(x-L)*scale)
term2s(x) = A*(2**(-D))*D*COS_T*exp(-(D+a)*L*scale)*cosh((D+a)*(x-L)*scale)
term2c(x) = B*(2**(-(2-D)))*(2-D)*COS_T*exp(-((2-D)+a)*L*scale)*cosh(((2-D)+a)*(x-L)*scale)
term3s(x) = A*(2**(-D))*D*(COS_T*COS_T*(D+1)-1)*exp(-(D+b)*L*scale)*cosh((D+b)*(x-L)*scale)
term3c(x) = B*(2**(-(2-D)))*(2-D)*(COS_T*COS_T*((2-D)+1)-1)*exp(-((2-D)+b)*L*scale)*cosh(((2-D)+b)*(x-L)*scale)
term4s(x) = 0#A*(2.0/3.0)*COS_T*D*(D+1)*(2*COS_T*COS_T*(D+2)-3)*exp(-(D+c)*L*scale)*cosh((D+c)*(x-L)*scale)
term4c(x) = 0#B*(2.0/3.0)*COS_T*(2-D)*((2-D)+1)*(2*COS_T*COS_T*((2-D)+2)-3)*exp(-((2-D)+c)*L*scale)*cosh(((2-D)+c)*(x-L)*scale)
term5(x) = 0#A*(1.0/6.0)*D*(D+1)*(4*COS_T*COS_T*COS_T*COS_T*(D*D+5*D+6)-12*COS_T*COS_T*(D+2)+3)*exp(-(D+f)*L*scale)*cosh((D+f)*(x-L)*scale)
term6(x) = 0#A*(1.0/15.0)*COS_T*D*(D+1)*(D+2)*(4*COS_T*COS_T*COS_T*COS_T*(D*D+7*D+12)-20*COS_T*COS_T*(D+3)+15)*exp(-(D+g)*L*scale)*cosh((D+g)*(x-L)*scale)

fit [2:L] term1s(x) + term1c(x) + term2s(x) + term2c(x) + term3s(x) + term3c(x) + term4s(x) + term4c(x) + term5(x) + term6(x) 'correlators.dat' using 1:2:3 via A,B,D

set label sprintf("Coefficients\nC_{0} = %e\nC_{1} = %e\nC_{2} = %e\nC_{3} = %e\n\nC_{/Symbol c}/C_{/Symbol s} = %f\n\n{/Symbol D} = %f\n\na = %f", 0.5*exp(-D*L*scale), D*COS_T*exp(-(D+a)*L*scale), D*(2*D+1)*COS_T*COS_T*exp(-(D+b)*L*scale), (2.0/3.0)*COS_T*D*(D+1)*(2*COS_T*COS_T*(D+2)-3)*exp(-(D+c)*L*scale), B/A, D, scale) at 2.0,0.9


   
plot [0:5] "correlators.dat" using 1:2:3 with yerrorbars t "Correlation function value", term1s(x) + term1c(x) + term2s(x) + term2c(x) + term3s(x) + term3c(x) + term4s(x) + term4c(x) + term5(x) t "sum", term1s(x) t "exp(-Dst)", term2s(x) t "exp(-t(Ds+1))",  term3s(x) t "exp(-t(Ds+2))", term4s(x) t "exp(-t(Ds+3))", term1c(x) t "exp(-Dct)", term2c(x) t "exp(-t(Dc+1))", term3c(x) t "exp(-t(Dc+2))", term4c(x) t "exp(-t(Dc+3))"#, term5(x) t "exp(-t(D+4))", term6(x) t "exp(-t(D+5))"
