#! /usr/bin/env bash

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
plot 'rb_prepare_data_solution.dat' u 1:a w l, 'rb_prepare_data_data.dat' u 1:a w p
#
set term wxt 1
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="K_ICF"
plot 'rb_prepare_data_solution.dat' u 1:a w l, 'rb_prepare_data_data.dat' u 1:a w p
#
set term wxt 2
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="s29"
s=sprintf("%s (%s)", a, "Glucose")
plot 'rb_prepare_data_solution.dat' u 1:a t s w l, 'rb_prepare_data_data.dat' u 1:a t s w p
#
set term wxt 3 
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="Insulin"
plot 'rb_prepare_data_solution.dat' u 1:a w l, 'rb_prepare_data_data.dat' u 1:a w p
#
set term wxt 4 
set grid
set xlabel "Time [h]"
set ylabel "Conc. [a.u.]"
a="K_urin"
plot 'rb_prepare_data_solution.dat' u 1:a w l, 'rb_prepare_data_data.dat' u 1:a w p
EOF
