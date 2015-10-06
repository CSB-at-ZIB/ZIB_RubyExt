#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/Model'
require_relative '../../lib/Experiment'
require_relative '../../lib/SysBioFit'

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

model = Model.new :reac_scheme, initvals    # t0, y0, par, plabel
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = $rtol
model.atol = $atol

# # puts "#{model.version}"
# n = 65
# tspan = (0..n).collect { |j| model.t0 + 150.0*j/n }
# # tspan = [ 0.0, 150.0 ]
#
# model.solve_ode tspan
# model.save_current_solution STDOUT
#
# exit -1


# n = 100
# tspan = (0..n).collect { |j| model.t0 + 150.0*j/n }
# # tspan = [ 0.0, 150.0 ]
# y0ode = model.y0ode
# pidx = [1,2,3]
# 
# model.solve_var tspan, [0.0, 0.0], model.par0, pidx
# model.save_current_solution STDOUT
# 
# 
# tp0, sol0 = model.solve_ode tspan, y0ode, model.par0
# 
# $jac = {}
# pidx.each_with_index do |ell,idx|
#    par = model.par0.clone
#    par[ell-1] += $eta 
#    STDERR.puts "par = #{par}"
#    tp, sol = model.solve_ode tspan, y0ode, par
#    tp0.each do |t|  
#       y = sol[t]
#       y0 = sol0[t]
#       $jac[t] ||= []
#       y0.each_with_index { |val,k| $jac[t] << (y[k] - y0[k])/$eta }
#    end
# end
# 
# ### puts "jac = #{$jac}"
# 
# ### STDOUT.printf "\n\n"
# printf "%-12s", "Timepoint"
# y0ode.each_with_index { |dmy,idx| printf "\t%12d", idx+1 }
# pidx.each_with_index do |ell,dmy| 
#   y0ode.each_with_index { |val,idx| printf "\t %05d/%05d", idx+1, ell }
# end
# printf "\n"
# tp0.each do |t|
#   printf "%.6e", t
#   y = sol0[t]
#   s = $jac[t]
#   y.each { |val| printf "\t%.6e", val }
#   s.each { |val| printf "\t%.6e", val }
#   printf "\n"
# end
# printf "\n\n"
# 
# exit -1

# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_NotIdentifiable_data.dat"

data = Experiment.new
data.load_csv fname

# ---------------------------------------------------------------------
# Parameter Estimation/Identification

pIniGuess = {
               "k1" => [ 0.6 , 0.0 ],
               "k2" => [ 0.2 , 0.0 ],
               "k3" => [ 1.5 , 0.0 ]

               # "k1" => [ 1.0,  0.0 ], 
               # "k2" => [ 0.5,  0.0 ],
               # "k3" => [ 0.08, 0.0 ] 
            }

pIniGuess = {
               "k1" => [ 5.569907e-01 , 0.0 ],
               "k2" => [ 1.723548e-01 , 0.0 ],
               "k3" => [ 1.391280e-01 , 0.0 ]
}

pIniGuess = {
               "k1" => [ 6.546714e-01 , 0.0 ],
               "k2" => [ 1.466385e-01 , 0.0 ],
               "k3" => [ 1.679909e-01 , 0.0 ]
}

nPar = 0
pIniGuess.keys.each do |key| 
  # next unless model.pId.index(key) 
  nPar += 1 if (model.pId.index(key) or model.y0Id.index(key)) 
end

mTotal = data.mdata.flatten.length
mFit   = mTotal - data.mweight.count(0.0)

fit = SysBioFit.new [nPar,mTotal,mFit]
fit.rtol = $ptol
fit.printlevel = 3
fit.pfname = "rb_Nlscon_with_NotIdentifiable_parameter.dat"
fit.sfname = "rb_Nlscon_with_NotIdentifiable_solution.dat"
fit.nitmodulo = 1
fit.nitmax = 35
fit.nonlin = 3
fit.jacgen = 1
fit.rwk = { "ajdel" => $eta, "cond" => 1e9 }

task = { model: model, data: data, guess: pIniGuess }

# fit.compute_sensitivity task

estim = {}
estim = fit.identify_par task

# ---------------------------------------------------------------------
# Result Output

labels = [].fill("n/a",0,nPar)
estim["pidx"].each_with_index do |j,idx|
  labels[idx] = model.pId[j-1] if j > 0
end

if estim.key?("par") then
  printf "\n               "
  # puts "                    k2           k1"
  labels.each { |str| printf "       %7s   ", str }
  printf "\n               "
  labels.each { |dmy| printf "  %14s ", "--------------" }
  printf "\n init. guess :"
  labels.each { |str| printf "  % 14.6f ", pIniGuess[str][0] }
  printf "\n               "
  labels.each { |dmy| printf "  %14s ", "--------------" }
  printf "\n confid_lower:"
  estim["rwk"]["xl"].each { |val| printf "  % 14.6f ", val }
  printf "\n  *p_estim*  :"
  estim["par"].each { |val| printf "  % 14.6f ", val }
  printf "\n confid_upper:"
  estim["rwk"]["xr"].each { |val| printf "  % 14.6f ", val }
  printf "\n\n"
  printf "incomp. kappa:     % .4e\n", estim["rwk"]["skap"]
  printf "achieved rtol:     % .4e\n", estim["rwk"]["prec"] 
  printf "\n"

  # puts "#{estim["iopt"]}"
  # puts "#{estim["rwk"]}"
end

puts "#{model.version}"
puts " "

