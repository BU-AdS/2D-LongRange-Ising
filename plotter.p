set autoscale    # scale axes automatically
set xtic auto    # set xtics automatically
set ytic auto    # set ytics automatically
#set format y ""
set xrange       [-50.00:10]
set yrange       [-50.00:10]
set output 'q=_Q__m2=_MASS__Cm2=_C_MASS__CN=_N_LATT_.eps'
set terminal postscript color enhanced
#set terminal dumb 120 40

set title "q=_Q_, Levels=_LEVELS_, _CENTRE_STRING_-centred, Dirichlet, m^{2} = _MASS_" 
set xlabel "log[2(1-|x|)(1-|y|)/((1-|x|)^{2}+(1-|y|)^{2}+|x-y|^{2})]" font "Times,16"
set ylabel "log[M^{-1}(i,j)]" font "Times,24"

deltaA = 1.0
NA = 1
deltaL = 1.0
NL = 1

actu(x) = deltaA*x + NA
latt(x) = deltaL*x + NL

cutoff = -100.0
upper  = -0.0

set print "crits_lev_LEVELS_.dat" append

fit [cutoff:upper] actu(x) 'lattice.dat' using (log($8)):(log($4)) via deltaA, NA
fit [cutoff:upper] latt(x) 'lattice.dat' using (log($8)):(log($2)) via deltaL, NL

print _MASS_,_MASS_/(deltaL*(deltaL-1)),exp(NA-NL),deltaL,deltaA

set label 1 sprintf("{/Symbol D}_{formula} = %.8f -> {/Symbol D}({/Symbol D}-1) = %.8f\n{/Symbol D}_{lattice} = %.8f -> {/Symbol D}({/Symbol D}-1) = %.8f\n\nC_{formula} = %.8f\nC_{lattice} = %.8f", deltaA, deltaA*(deltaA-1), deltaL, deltaL*(deltaL-1), NA, NL) at -15,(-30) font "courier,14"

set label 2 sprintf("Scaling Factors\n{/Symbol g}_{m^{2}} = %.8f  {/Symbol g}_{N} = %.8f", _C_MASS_, _N_LATT_) at -15,(-40) font "courier,14"

set label 3 sprintf("Try %.12f %.8f", deltaA*(deltaA-1)/(deltaL*(deltaL-1)), exp(NA-NL)) at -15,(-45) font "courier,14"

plot "lattice.dat" using (log($8)):(log($4)) t 'formula',  "lattice.dat" using (log($8)):(log($2)) t 'lattice', latt(x) t '', actu(x) t ''
