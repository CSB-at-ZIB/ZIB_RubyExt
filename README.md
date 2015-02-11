# ZIB_RubyExt
## Wrapping ancient FORTRAN routines to Ruby via Ruby Extensions

The package ZIB_RubyExt offers a really simple and down-to-earth approach
to make use of ancient FORTRAN routines in Ruby.  Currently, an error-oriented,
globalised Gauss-Newton algorithm (with possible nonlinear constraints 
treatment), and a linearly-implicit ODE solver with extrapolation, NLSCON 
and LIMEX, respectively, are provided as FORTRAN examples.  These Ruby 
Extension classes are used as base classes for derived Ruby classes to 
solve parameter estimation/identification problems in ODE systems, as they  
appear often in systems biology.


## License

This software package is released under GPL 2.0, see LICENSE file.


## Requirements

* Ruby (1.9.3p0 or higher patch levels)
* (any) compatible Ruby-Dev Package
* gfortran (>= 4.8.2)
* g++ (>= 4.8.2)

The package has been tested and runs with Ruby 1.9.3p0 (or higher 
patch levels).  For the make process to finish successully, any kind 
of ruby-dev package is needed (e.g. ruby-dev_1.9.3.4_all.deb for 
Ubuntu 14.04.2 LTS; for other Linux distributions probably similar).

The makefiles expect that gfortran is available, as well as g++ 
(both >= 4.8.2).

### Optional

The scripts displaying all the resulting tables of tons of numbers use
* [gnuplot](http://www.gnuplot.info)



## Installation

````
[.../path_to/ZIB_RubyExt]$ make clean; make
````

NOTE: The first time this command is invoked could take slighly a bit longer
since a libsbml is build.  For the sake of completeness, and ease of use, 
this (optional) step has been included as it is needed for the compiler 
``sbml2fortran´´.


