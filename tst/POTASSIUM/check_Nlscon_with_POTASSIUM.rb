#! /usr/bin/env ruby

# =====================================================================

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/ModelDL'
require_relative '../../lib/Experiment'
require_relative '../../lib/SysBioFit'

# =====================================================================



# ---------------------------------------------------------------------
# Model/dynamic-load ODE

model = ModelDL.new 
model.t0 = 5.0
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = 1.0e-9
model.atol = 1.0e-9
# model.monitor = 1
# puts "#{model.inspect}"
# puts "#{model.version}" # see at the end of this script



# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_prepare_data_data.dat"

data = Experiment.new 
data.load_data fname

#
# ... check differences between experimental data and model simulation
#
midx = []
data.mdata.each do |ary|
  row = []
  ary.each_with_index do |idx,j|
    next if idx.is_a?(String)
    label = data.mlabel[j]
    k = model.yId.index(label)
    row << k
  end
  midx << row
end

tspan = [model.t0].concat(data.mtime.sort).uniq
ts, solution = model.solve_ode(tspan)

cnt = 0
printf("\nExperimental data _and_ simulated data : \n")
data.mtime.length.times do |k|
  y = solution[ data.mtime[k] ]
  midx[k].each_with_index do |idx,j|
    val = data.mdata[k][j]
    printf("  %4d:   %12.7f (SD: %.3e)     %12.7f     diff: %.6e\n", 
             cnt, val, data.mweight[k][j], y[idx], (val - y[idx]).abs)
    cnt += 1
  end
end
printf("\n\n\n")
#
# end checking!
#



# ---------------------------------------------------------------------
# Parameter Estimation/Identification

pIniGuess = {
                "global_p2"   =>  [ 34.0   ,  0.0 ],   #  par(2)
                "global_p4"   =>  [  0.0002,  0.0 ],   #  par(4)
                "global_p5"   =>  [ 14.04  ,  0.0 ],   #  par(5)
                "global_p6"   =>  [  0.05  ,  0.0 ],   #  par(6)
                "global_p8"   =>  [ 14.91  ,  0.0 ],   #  par(8)
                "global_p9"   =>  [  9.99  ,  0.0 ],   #  par(9)
                "global_p13"  =>  [  6.01  ,  0.0 ],   #  par(13)
                "global_p16"  =>  [  3.0   ,  0.0 ],   #  par(16)
                "global_p18"  =>  [  0.08  ,  0.0 ],   #  par(18)
                "global_p19"  =>  [  0.1373,  0.0 ],   #  par(19)
                "global_p21"  =>  [  0.1085,  0.0 ],   #  par(21)
                "global_p25"  =>  [  0.2898,  0.0 ],   #  par(25)
                "global_p26"  =>  [  2.15  ,  0.0 ],   #  par(26)
 
                "global_p44"  =>  [  1.2   ,  0.0 ],   #  par(42)
                "global_p45"  =>  [  3.0   ,  0.0 ],   #  par(43)
                "global_p46"  =>  [  0.0854,  0.0 ],   #  par(44)
                "global_p47"  =>  [  1.0   ,  0.0 ],   #  par(45)
                "global_p48"  =>  [  0.015 ,  0.0 ],   #  par(46)
                "global_p49"  =>  [  0.985 ,  0.0 ],   #  par(47)
                "global_p50"  =>  [  7.0   ,  0.0 ]    #  par(48)
            }


nPar = 0
pIniGuess.keys.each do |key|
  # next unless model.pId.index(key)
  nPar += 1 if (model.pId.index(key) or model.y0Id.index(key))
end
# nPar   = pidx.length
mTotal = data.mdata.flatten.length - data.mdata.flatten.count("n/a")
mFit   = mTotal - data.mweight.count(0.0)


fit = SysBioFit.new [nPar,mTotal,mFit]
fit.rtol = 1.0e-2
fit.printlevel = 3
fit.pfname = "rb_Nlscon_with_POTASSIUM_parameter.dat"
fit.sfname = "rb_Nlscon_with_POTASSIUM_solution.dat"
#fit.nitmodulo = 5 # every nitmodulo-th solution is written to sfname 
fit.nitmax = 30
fit.nonlin = 4
fit.rwk = { "cond" => 5.0e+5 } 


id_task = { model: model, data: data,  guess: pIniGuess }
           # pidx: pidx, guess: guess, pscal: pscal }

         
estim = fit.identify_par id_task 



# ---------------------------------------------------------------------
# Result Output

puts " "
puts "incomp. kappa:  #{'% .4e' % estim["rwk"]["skap"]}"
puts "achieved rtol:  #{'% .4e' % estim["rwk"]["prec"]}" 
puts " "

# puts "#{estim["iopt"]}"
# puts "#{estim["rwk"]}"

# puts "#{estim["t"]}"
# puts "#{estim["y"]}"

puts "       pidx: #{fit.pidx}"
puts "   pinitial: #{fit.set_current_par(fit.pguess)}"
puts "  model par: #{model.par0}"
puts "  model pId: #{model.pId}"
puts "  final par: #{fit.set_current_par(fit.x)}"

puts " "
puts "model version: #{model.version}"
puts " "

