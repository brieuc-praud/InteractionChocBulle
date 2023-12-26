set fit quiet
set fit logfile "/dev/null"
set key right bottom

ERRFILE = "error.dat"

#fit f(x) ERRFILE u (log($1)):(log($3)) via b

set logscale xy
set format x "%2.0t.10^{%L}"
set format y "%2.0t.10^{%L}"

set title "Error"
plot ERRFILE u 1:2 w lp t "density"
replot ERRFILE u 1:3 w lp t "velocity_x"
replot ERRFILE u 1:4 w lp t "velocity_y"
replot ERRFILE u 1:5 w lp t "energy"
f(x)=1.5*GPVAL_DATA_Y_MAX/GPVAL_DATA_X_MIN * x
replot f(x) t "slope 1"
