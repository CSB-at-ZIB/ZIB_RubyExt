#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/ModelDL'

# ---------------------------------------------------------------------
# Command line

$rtol = 1.0e-6
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
  puts "          * rtol = atol = 1.0e-6"
  puts " "
  exit -1
end
$rtol = 1.0e-6 if $rtol <= 0.0
$atol = 1.0e-6 if $atol <= 0.0

model = ModelDL.new 
model.t0 = 0.0
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = $rtol
model.atol = $atol
# model.monitor = 1

puts ""
puts "ModelDL settings:"
puts "pId  : #{model.pId}"
puts "par0 : #{model.par0}"
# puts "#{model.version}"

# ---------------------------------------------------------------------
# Parameter Estimation/Identification


pIniGuess = {
   "global_p_019_001"  =>  [    0.958 ,  1.0 ],  # GynCycle Param 61
   "global_p_020_001"  =>  [    0.925 ,  1.0 ],  # GynCycle Param 62
   "global_p_021_001"  =>  [    0.7576,  1.0 ],  # GynCycle Param 63
   "global_p_022_001"  =>  [    0.61  ,  1.0 ],  # GynCycle Param 64
   "global_p_023_001"  =>  [    0.543 ,  1.0 ]   # GynCycle Param 65

#  "global_p_024_001"  =>  [   51.558 ,  1.0 ],  # GynCycle Param 69

#  "global_p_024_002"  =>  [    2.0945,  1.0 ],  # GynCycle Param 70
#  "global_p_024_003"  =>  [    9.28  ,  1.0 ],  # GynCycle Param 71
#  "global_p_024_004"  =>  [ 6960.53  ,  1.0 ],  # GynCycle Param 72
#  "global_p_024_005"  =>  [    0.972 ,  1.0 ],  # GynCycle Param 73
#  "global_p_024_006"  =>  [ 1713.71  ,  1.0 ],  # GynCycle Param 74
#  "global_p_024_007"  =>  [ 8675.14  ,  1.0 ],  # GynCycle Param 75
#  "global_p_024_008"  =>  [    5.235 ,  1.0 ]   # GynCycle Param 76
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

tspan = [0.0,30.0]
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

puts ""
puts "ModelDL: model.solve_var called with:"
puts "pidx : #{pidx}"
puts "par0 : #{model.par}"

model.solve_var tspan, y0, par, pidx

# ---------------------------------------------------------------------
# Result Output

sout = File.open("rb_ModelDL_vareq_PAEON_V2_solution.dat","w")
model.save_current_solution sout
sout.close

puts " "
puts "#{model.version}"
puts " "

