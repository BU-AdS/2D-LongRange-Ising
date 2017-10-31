set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [-0.5:6.5]
#set yrange       [-0:2.5]
set output 'test.eps'
set terminal postscript color enhanced
#set terminal dumb 160 40

set title "q=7, L=4, T=20"
set xlabel "i" font "Times-Roman,24"
set ylabel "eval" font "Times-Roman,24"

f(x) = A + B*sqrt(x)
#g(x) = a + b*x**n

fit [0:2500] f(x) 'test.dat' using 1:2 via A,B

set label 1 sprintf("f(x) = A + B*sqrt(x) from x=0..2500\n\nA=%f, B=%f", A, B) at 1500,(6.0) font "courier,14"

plot [1:4639] "test.dat" using 1:2 t "spectrum", f(x) t "fit" 
