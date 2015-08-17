#! /usr/bin/env bash

tp="0 -1"
if [ $# -gt 0 ]; then
  tp=${1}
fi

fnsol="rb_Nlscon_with_POTASSIUM_solution.dat"
fndat="rb_prepare_data_data.dat"
#fndat="rb_Nlscon_with_POTASSIUM_data.dat"

#
# Use with care! *All* plots are closed by this command:
#
#killall /usr/bin/gnuplot 2>/dev/null

/usr/bin/gnuplot -p << EOF
set term wxt 0
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="K_ECF"
plot for [k in "${tp}"] '${fnsol}' i sprintf("== Iter %4s ==", k) u 1:a t sprintf("Iter %4s",k) w l lw 3, '${fndat}' u 1:a t col w p
#
set term wxt 1
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="K_ICF"
plot for [k in "${tp}"] '${fnsol}' i sprintf("== Iter %4s ==", k) u 1:a t sprintf("Iter %4s",k) w l lw 3, '${fndat}' u 1:a t col w p
#
set term wxt 2
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="s29"
s=sprintf("%s (%s)", a, "Glucose")
plot for [k in "${tp}"] '${fnsol}' i sprintf("== Iter %4s ==", k) u 1:a t sprintf("Iter %4s",k) w l lw 3, '${fndat}' u 1:a t s w p
#
set term wxt 3 
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="Insulin"
plot for [k in "${tp}"] '${fnsol}' i sprintf("== Iter %4s ==", k) u 1:a t sprintf("Iter %4s",k) w l lw 3, '${fndat}' u 1:a t col w p
#
set term wxt 4 
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="K_urin"
plot for [k in "${tp}"] '${fnsol}' i sprintf("== Iter %4s ==", k) u 1:a t sprintf("Iter %4s",k) w l lw 3, '${fndat}' u 1:a t col w p
EOF
