#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/Model'
require_relative '../../lib/Experiment'
# require_relative '../../lib/SysBioFit'

# ---------------------------------------------------------------------
# Command line

$rtol = 1.0e-8
$atol = 1.0e-9

$ptol = 1.0e-1
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
  puts "          * ptol = 1.0e-1"
  puts "          * eta = sqrt(10.0*rtol) ; rtol = #{'%.1e' % $rtol}"
  puts " "
  exit -1
end
$ptol = 1.0e-1 if $ptol <= 0.0
$eta = Math::sqrt(10.0*$rtol) if $eta <= 0.0

# ---------------------------------------------------------------------
# ODE system

def control(t)
  (1.0 + 0.75*Math::sin(0.37*Math::PI*t))*Math::exp(-0.025*(t-20.0).abs**1.25)
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

#

def linspace(a,b,num=100)
  num -= 1
  arr = [*0..num]
  arr.collect! {|j| a + j.to_f*(b-a)/num}
end

# ---------------------------------------------------------------------
# Model/ODE setup

initvals = {
                   t0:  0.0 , 
                   y0: [ 0.0  , 0.0  ] ,
              y0label: [ "A0" , "B0" ] ,

                  par: [  1.2 ,  0.08 ,  0.5  ] ,
               plabel: [  "k1",   "k2",  "k3" ] ,

                  jac: :jacobian
           }

$model = Model.new :reac_scheme, initvals    # t0, y0, par, plabel
$model.hmax = 0.0
$model.inistep = 1.0e-4
$model.rtol = $rtol
$model.atol = $atol

# puts "#{model.version}"
# exit -1

# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_NotIdentifiable_data.dat"

data = Experiment.new
data.load_csv fname

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

tspan = [$model.t0].concat(data.mtime.sort).uniq
meas = data.mdata.flatten
par = $model.par0.clone

llpar = par.collect { |x| Math::log(x/10.0) }
lupar = par.collect { |y| Math::log(10.0*y) }

pvec = []
par.length.times { |k| pvec[k] = linspace(llpar[k],lupar[k],20) }

# pvec.each { |arr| printf "\t % .6e\t % .6e\n", Math::exp(arr[0]), Math::exp(arr[-1]) }
#
# exit -1

printf "Distance"
par.each_with_index {|x,idx| printf "\t   %5d", idx+1}
printf "\n"
pvec[0].each do |x|
  pvec[1].each do |y|
    pvec[2].each do |z|
       par = [Math::exp(x),Math::exp(y),Math::exp(z)]
       t, sol = $model.solve_ode tspan, $model.y0, par
       l2sum = 0.0
       meas.length.times do |j|
          tj = tspan[j+1]
          yj = sol[tj][1] if sol.has_key?(tj)
          l2sum += (yj - meas[j])**2
       end
       printf "% .6e", Math::sqrt(l2sum)/meas.length
       par.each {|val| printf "\t % .6e", val}
       printf "\n"
    end
  end
end


