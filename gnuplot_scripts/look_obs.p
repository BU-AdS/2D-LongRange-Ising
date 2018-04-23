set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
#set xrange       [0:10]
#set yrange       [-2:0]
set terminal pdf
set output "plot_obs.pdf"
#set terminal dumb 120 40

plot [0:1000] "observables.dat" using 1:4 t "" ps 1 pt 0
