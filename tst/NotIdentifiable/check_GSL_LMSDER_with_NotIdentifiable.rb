#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require("gsl")
include GSL::MultiFit

require_relative '../../lib/Model'
require_relative '../../lib/Experiment'
require_relative '../../lib/SysBioFit'

# ---------------------------------------------------------------------
# Command line

$rtol = 1.0e-8
$atol = 1.0e-9

$ptol = 1.0e-2
$eta = Math::sqrt(10.0*$rtol)

if ARGV.length > 0 then
  $ptol = Float(ARGV[0]) rescue 0.0
  if ARGV.length > 1 then
    $eta = Float(ARGV[1]) rescue 0.0
  end
else
  puts " "
  puts "Usage: ./#{File.basename(__FILE__)} ptol [eta]"
  puts " "
  puts "         Default values"
  puts "         --------------"
  puts "          * ptol = 1.0e-2"
  puts "          * eta = sqrt(10.0*rtol) ; rtol = #{'%.1e' % $rtol}"
  puts " "
  exit -1
end
$ptol = 1.0e-2 if $ptol <= 0.0
$eta = Math::sqrt(10.0*$rtol) if $eta <= 0.0

# ---------------------------------------------------------------------
# ODE system

def control(t)
  (1.0 + 0.75*Math::sin(0.37*Math::PI*t))*
            Math::exp(-0.025*(t-20.0).abs**1.25)
end

def reac_scheme(t,y,par)

  k1, k2, k3 = par

  r0 = k1 * y[0]
  r1 = k2 * y[0]
  r2 = k3 * y[1]

  u = control(t)
  
  [ 
    - r0 - r1 + r2 + u,
      r0 - r2 
  ]

end

#

def jacobian(t,y,par,pidx)

  dy = reac_scheme(t,y,par) # y is in fact too long, but this does not matter! 

  nspe = dy.length
  npar = par.length
  k1, k2, k3 = par

  drdy = [ [  k1 , 0.0 ] ,
           [  k2 , 0.0 ] ,
           [ 0.0 ,  k3 ] ]

  drdp = [ [ y[0] ,  0.0 ,  0.0 ] ,
           [  0.0 , y[0] ,  0.0 ] ,
           [  0.0 ,  0.0 , y[1] ] ]

  s = y[nspe..-1]  # m = 2 species (y[0], y[1]), and q = 3 parameters
                   #  ==>  s ~ y[m = 2], ..., y[(m+m*q-1) = (2 + 2*3 - 1) = 7]

  fy = [ [ - drdy[0][0] - drdy[1][0] + drdy[2][0] , - drdy[0][1] - drdy[1][1] + drdy[2][1] ] ,
         [   drdy[0][0] - drdy[2][0]              ,   drdy[0][1] - drdy[2][1]              ] ]

  fp = [ [ - drdp[0][0] - drdp[1][0] + drdp[2][0] , - drdp[0][1] - drdp[1][1] + drdp[2][1] , - drdp[0][2] - drdp[1][2] + drdp[2][2] ] ,
         [   drdp[0][0] - drdp[2][0]              ,   drdp[0][1] - drdp[2][1]              ,   drdp[0][2] - drdp[2][2]              ] ]

  dS = []
  pidx.each_with_index do |ell,idx|
    # next if ell < 0
    nspe.times do |j|
      sum = 0.0
      sum = fp[j][ell-1] if ell > 0
      nspe.times { |nu| sum += fy[j][nu]*s[nu + nspe*idx] }
      dS[j + nspe*idx] = sum
    end
  end

  [ dy , dS ].flatten!

end

# ---------------------------------------------------------------------
# Model/ODE setup

initvals = {
                   t0:  0.0 , 
                   y0: [ 0.0  , 0.0   ] ,
              y0label: [ "A0" , "B0"  ] ,

                  par: [  1.2 ,  0.08 ,  0.5  ] ,
               plabel: [  "k1",   "k2",  "k3" ] ,

                  jac: :jacobian
           }

$model = Model.new :reac_scheme, initvals    # t0, y0, par, plabel
$model.hmax = 0.0
$model.inistep = 1.0e-4
$model.rtol = $rtol
$model.atol = $atol

# # puts "#{$model.version}"
@pidx = []
(1..$model.par0.length).each {|j| @pidx << j }

# ---------------------------------------------------------------------

# x: Vector, list of the parameters to determine
# t, y, sigma: Vectors, observational data
# f: Vector, function to minimize
procf = Proc.new { |x, t, y, sigma, f|
  par = x.to_a
  n = t.size
  tspan = [$model.t0].concat(t.to_a.sort).uniq
  tp, sol = $model.solve_ode tspan, $model.y0ode, par
  for j in 0...n do
    tj = t[j]
    yj = sol[tj][1]
    f[j] = (yj - y[j])/sigma[j]
  end
}

# jac: Matrix, Jacobian
procdf = Proc.new { |x, t, y, sigma, jac|
  par = x.to_a
  m = $model.y0ode.size
  n = t.size
  tspan = [$model.t0].concat(t.to_a.sort).uniq
  tp, sens = $model.solve_var tspan, $model.y0ode, par, @pidx
  for j in 0...n do
    tj = t[j]
    sj = sigma[j]
    @pidx.each do |k|
      jac.set(j, k-1, sens[tj][m+(m*(k-1)+1)] / sj )
    end
  end
}

f = Function_fdf.alloc(procf, procdf, 2)

# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_NotIdentifiable_data.dat"

data = Experiment.new
data.load_csv fname
n = data.mdata.flatten.length

# ---------------------------------------------------------------------

t = GSL::Vector.alloc(n)
y = GSL::Vector.alloc(n)
sigma = GSL::Vector.alloc(n)
for j in 0...n do
  t[j] = data.mtime[j].to_f
  y[j] = data.mdata.flatten[j]
  sigma[j] = data.mweight.flatten[j]
end

f.set_data(t, y, sigma)
#x = GSL::Vector.alloc(0.6, 0.2, 1.5)    # initial guess
# x = GSL::Vector.alloc(1.0, 0.5, 0.08)    # initial guess
# x = GSL::Vector.alloc(5.569907e-01, 1.723548e-01, 1.391280e-01)    # initial guess
x = GSL::Vector.alloc(6.546714e-01 , 1.466385e-01, 1.679909e-01)    # initial guess
np = 3

solver = FdfSolver.alloc(FdfSolver::LMSDER, n, np)

solver.set(f, x)

printf("n = %d (data size)\n", n)
printf "inichi/n = %e\n", solver.f.dnrm2/n

iter = 0
solver.print_state(iter)
begin
  iter += 1
  status = solver.iterate
  solver.print_state(iter)
  status = solver.test_delta($ptol, $ptol)
end while status == GSL::CONTINUE and iter < 500

covar = solver.covar(0.0)
position = solver.position
chi = solver.f.dnrm2
dof = n - np
printf("dof = %d (= %d - %d )\n", dof, n, np)
printf("chi2/dof = %e\n", chi**2/dof)
printf("k1 = %.5f +/- %.5f\n", position[0], chi*Math::sqrt(covar[0,0]/dof))
printf("k2 = %.5f +/- %.5f\n", position[1], chi*Math::sqrt(covar[1,1]/dof))
printf("k3 = %.5f +/- %.5f\n", position[2], chi*Math::sqrt(covar[2,2]/dof))

tspan = [ $model.t0, t.to_a[-1] ]

$model.solve_ode tspan, $model.y0ode, position.to_a 
$model.save_current_solution STDERR

