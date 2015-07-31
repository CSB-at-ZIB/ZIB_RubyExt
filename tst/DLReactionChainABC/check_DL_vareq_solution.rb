#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/ModelDL'

# ---------------------------------------------------------------------
# Command line

$rtol = 1.0e-8
$atol = $rtol

if ARGV.length > 0 then
  $rtol = Float(ARGV[0]) rescue 0.0
  if ARGV.length > 1 then
    $atol = Float(ARGV[1]) rescue 0.0
  end
else
  puts " "
  puts "Usage: ./#{File.basename(__FILE__)} rtol [atol]"
  puts " "
  puts "         Default values"
  puts "         --------------"
  puts "          * rtol = atol = 1.0e-8"
  puts " "
  exit -1
end
$rtol = 1.0e-8 if $rtol <= 0.0
$atol = 1.0e-8 if $atol <= 0.0

model = ModelDL.new 
model.t0 = 0.0
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = $rtol
model.atol = $atol
# model.monitor = 1

# puts "#{model.version}"

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

tspan = [0.0,5.0]
y0 = model.y0ode
par = model.par0
pidx = []
pIniGuess.keys.each do |key|
   idx = model.pId.index(key) 
   if idx then
      pidx << idx+1
      par[ idx ] = pIniGuess[key][0]
   end
end

model.solve_var tspan, y0, par, pidx

# ---------------------------------------------------------------------
# Result Output

sout = File.open("rb_DL_vareq_solution.dat","w")
model.save_current_solution sout
sout.close

puts "#{model.version}"
puts " "

