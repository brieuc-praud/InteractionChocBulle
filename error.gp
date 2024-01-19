set fit quiet
set fit logfile "/dev/null"
set key right bottom

ERRFILE = "error.dat"

set logscale xy
set format x "%2.0t.10^{%L}"
set format y "%2.0t.10^{%L}"

set title "Error"

f(x)=a*x+b
fit f(x) ERRFILE u (log($1)):(log($2)) via a,b
plot exp(f(log(x))) lc rgb 'gray' dt 4 t sprintf("slope %.2f", a)

replot ERRFILE u 1:2 w lp t "density"
replot ERRFILE u 1:3 w lp t "velocity_x"
replot ERRFILE u 1:4 w lp t "velocity_y"
replot ERRFILE u 1:5 w lp t "energy"

set xrange [0.9*GPVAL_X_MIN:1.1*GPVAL_X_MAX]
set yrange [1.5*GPVAL_Y_MIN:1.5*GPVAL_Y_MAX]
replot
