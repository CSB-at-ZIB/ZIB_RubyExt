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
model.t0 = -42.5
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = 1.0e-9
model.atol = 1.0e-9
# model.monitor = 1
# puts "#{model.inspect}"
puts "#{model.version}" # see at the end of this script



# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_PAEON_V2_data.dat"

data = Experiment.new 
data.load_data fname

#
# ... check differences between experimental data and model simulation
#
fobs  = data.mdata.flatten.delete_if { |val| val.is_a?(String) }
fscal = data.mweight.flatten.delete_if { |val| val.is_a?(String) }
midx  = []
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
  next if midx[k] == []
  y = solution[ data.mtime[k] ]
  midx[k].each_with_index do |idx,j|
    val = fobs[cnt]
    printf("  %4d:   %12.7f (SD: %.3e)     %12.7f     diff: %.6e\n", 
             cnt, val, fscal[cnt], y[idx], (val - y[idx]).abs)
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

# ---------------------------------------------------------------------
nPar = 0
pIniGuess.keys.each do |key|
  # next unless model.pId.index(key)
  nPar += 1 if (model.pId.index(key) or model.y0Id.index(key))
end
# nPar   = pidx.length
mTotal = data.mdata.flatten.length - data.mdata.flatten.count("n/a")
mFit   = mTotal - data.mweight.count(0.0)
# ---------------------------------------------------------------------

nlscon = SysBioFit.new [nPar,mTotal,mFit]
nlscon.rtol = 1.0e-2
nlscon.printlevel = 3
nlscon.pfname = "rb_Nlscon_with_PAEON_V2_parameter.dat"
nlscon.sfname = "rb_Nlscon_with_PAEON_V2_solution.dat"
#nlscon.nitmodulo = 5 # every nitmodulo-th solution is written to sfname 
nlscon.nitmax = 45
nlscon.nonlin = 3
nlscon.rwk = { "cond" => 1.0e+9 } 


current_task = { model: model, data: data,  guess: pIniGuess }
                # pidx: pidx, guess: guess, pscal: pscal }

estim = {}         
estim = nlscon.identify_par current_task 

# ---------------------------------------------------------------------
# Result Output

if estim.has_key?("par") then
  printf "\n"
  printf "   %22s     %12s   %12s     %12s\n", 
         "Parameter", "confid_lower", "*p_estim*", "confid_upper"
  printf "   %22s     %12s   %12s     %12s\n", 
         "---------", "------------", "---------", "------------"

  estim["pidx"].each_with_index do |j,idx|
    label = (j > 0) ? model.pId[j-1] : model.y0Id[-j-1] 
    printf "  %22s     % 12.6f    % 12.6f    % 12.6f\n", 
            label[0..21], 
            estim["rwk"]["xl"][idx], 
            estim["par"][idx], 
            estim["rwk"]["xr"][idx]
  end
  printf "\n"

  puts " "
  puts "incomp. kappa:  #{'% .4e' % estim["rwk"]["skap"]}"
  puts "achieved rtol:  #{'% .4e' % estim["rwk"]["prec"]}"
  puts " "
end

# puts "#{estim["iopt"]}"
# puts "#{estim["rwk"]}"

# puts "#{estim["t"]}"
# puts "#{estim["y"]}"

# puts "       pidx: #{nlscon.pidx}"
# puts "   pinitial: #{nlscon.set_current_par(nlscon.pguess)}"
# puts "  model par: #{model.par0}"
# puts "  model pId: #{model.pId}"
# puts "  final par: #{nlscon.set_current_par(nlscon.x)}"

puts " "
puts "#{model.version}"
puts " "

# exit

# ---------------------------------------------------------------------
# Sensitivity Matrix Output (w.r.t. initial guess)

sens = nlscon.compute_sensitivity current_task, "l2"


# puts "sens['pidx'] = #{sens["pidx"]}"
puts " "
puts "                ------------------------------------------"
puts "                  sens['jcnrm'] (each pidx normed by #{sens["nrm"]})  "
puts "                ------------------------------------------"
puts " "

channel = [ $stdout , File.open("rb_sensitivity.dat", "w") ]
channel.each do |out|
  sens["jacobian"].each_with_index do |ary,j| 
    out.printf("# pidx = % 3d :  l2 = %e\n", 
                 sens["pidx"][j], sens["jcnrm"][j])
    ary.each_with_index do |val,n| 
       out.printf("  % .3e , ", val) 
       out.printf("\n") if (n+1).modulo(5) == 0
    end
    out.printf("\n\n")
  end
  out.printf("\n")
  out.printf("# midx : \"(k, j)\"  [indices of the measurement pair (t_k, y_j)]\n")
  cnt = 0
  data.mtime.length.times do |k|
    next if midx[k]==[]
    midx[k].each_with_index do |idx,j|
      cnt += 1
      out.printf(" (%4d, %3d) , ", k, idx)
      out.printf("\n") if cnt.modulo(5) == 0
    end
  end
  out.printf("\n\n")
  out.printf("# %s\n", "#{model.version}")
  out.printf("\n")
end

channel.shift
channel.each { |out| out.close }
