set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
#set yrange       [0:15]
set terminal pdf
set output "plot_crits.pdf"
#set terminal dumb 120 40

SHscale=20
Bscale=100

plot [2.0:2.5] "test.dat" using (1.0/($3)):4:5 with errorbars t "Susceptibiliy" ps 1 pt 1, "test.dat" using (1.0/($3)):(SHscale*($6)):(SHscale*($7)) with errorbars t "Specific Heat (Scaled)" ps 1 pt 2 , "test.dat" using (1.0/($3)):(Bscale*($8)):(Bscale*($9)) with errorbars t "Binder Cumulant (Scaled)" ps 1 pt 3
