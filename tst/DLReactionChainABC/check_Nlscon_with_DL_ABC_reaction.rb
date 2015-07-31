#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/ModelDL'
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
$eta = Math::sqrt(10.0*$rtol) if $ptol <= 0.0

# ---------------------------------------------------------------------
# ODE system

def abc(t,y,par)

  k1,k2 = par

  r1 = k1 * y[0]
  r2 = k2 * y[1]

  # dy0     dy1      dy2
  [ -r1 , +r1 - r2 , +r2 ]

end

# ---------------------------------------------------------------------
# Model/ODE setup

model = ModelDL.new
model.t0 = 0.0
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = $rtol
model.atol = $atol

puts "#{model.version}"
### puts "#{model.pId}"

# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_DL_ABC_data.dat"

data = Experiment.new 
data.load_csv fname

# ---------------------------------------------------------------------
# Parameter Estimation/Identification

pIniGuess = {
               "global_k2" => [ 3.5 , 0.0 ],
               "global_k1" => [ 0.3 , 0.0 ]
            }

nPar = 0
pIniGuess.keys.each do |key| 
  # next unless model.pId.index(key) 
  nPar += 1 if (model.pId.index(key) or model.y0Id.index(key)) 
end
# pidx   = [   2,   1 ]
# guess  = [ 3.5, 0.3 ]
# pscal  = [ 1.0, 1.0 ]

# nPar   = pidx.length
mTotal = data.mdata.flatten.length
mFit   = mTotal - data.mweight.count(0.0)

fit = SysBioFit.new [nPar,mTotal,mFit]
fit.rtol = $ptol
fit.printlevel = 3
fit.pfname = "rb_Nlscon_with_DL_ABC_parameter.dat"
fit.sfname = "rb_Nlscon_with_DL_ABC_solution.dat"
fit.nitmodulo = 1
fit.nitmax = 195
fit.nonlin = 3
fit.jacgen = 1
fit.rwk = { "ajdel" => $eta } 
# fit.rwk = { "cond" => 1.0e-3 } 

id_task = { model: model, data: data,  guess: pIniGuess }
           # pidx: pidx, guess: guess, pscal: pscal }
         
estim = fit.identify_par id_task 

### sens = fit.compute_sensitivity id_task

# ---------------------------------------------------------------------
# Result Output

labels = [].fill("n/a",0,nPar)
estim["pidx"].each_with_index do |j,idx| 
  labels[idx] = model.pId[j-1] if j > 0
end

puts " "
# puts "                    k2           k1"
puts "                  #{'%4s         %4s' % labels}"
puts "                ----------   ----------"
puts " confid_lower:  #{'% .6f    % .6f' % estim["rwk"]["xl"]}"
puts "  *p_estim*  :  #{'% .6f    % .6f' % estim["par"]}"
puts " confid_upper:  #{'% .6f    % .6f' % estim["rwk"]["xr"]}"
puts " "
puts "incomp. kappa:       #{'% .4e' % estim["rwk"]["skap"]}"
puts "achieved rtol:       #{'% .4e' % estim["rwk"]["prec"]}"
puts " "

# puts "#{estim["iopt"]}"
# puts "#{estim["rwk"]}"

# puts "#{estim["t"]}"
# puts "#{estim["y"]}"

puts "#{model.version}"
puts " "

