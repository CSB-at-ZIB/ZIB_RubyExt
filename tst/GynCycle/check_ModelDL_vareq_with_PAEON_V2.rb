#! /usr/bin/env ruby

begin
  require 'gsl'
  include GSL
  $gsl_avail = true
rescue LoadError
  puts " "
  puts "Hint: Unable to load ruby's gsl library."
  puts "      Probably try to install it by '[sudo] gem install gsl'"
  puts "      to get this script fully functional. Some parts"
  puts "      of it are switched off by now."
  puts " "
  $gsl_avail = false
end

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/ModelDL'

# ---------------------------------------------------------------------
# Command line

$rtol = 1.0e-6
$atol = $rtol
$delp = 0.25

if ARGV.length > 0 then
  $rtol = Float(ARGV[0]) rescue 0.0
  if ARGV.length > 1 then
    $atol = Float(ARGV[1]) rescue 0.0
    if ARGV.length > 2 then 
      $delp = Float(ARGV[2]) rescue 0.0
    end
  end
else
  puts " "
  puts "Usage: ./#{File.basename(__FILE__)} rtol [atol [dp]]"
  puts " "
  puts "         Default values"
  puts "         --------------"
  puts "          * rtol = atol = 1.0e-6"
  puts "          * dp = 0.25 (percent of p)"
  puts " "
  exit -1
end
$rtol = 1.0e-6 if $rtol <= 0.0
$atol = 1.0e-6 if $atol <= 0.0
$delp = 0.25 if $delp <= 0.0

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
   "global_p_021_001"  =>  [    0.7567,  1.0 ],  # GynCycle Param 63
   "global_p_022_001"  =>  [    0.61  ,  1.0 ],  # GynCycle Param 64
   "global_p_023_001"  =>  [    0.543 ,  1.0 ],  # GynCycle Param 65

#  "global_p_024_001"  =>  [   51.558 ,  1.0 ],  # GynCycle Param 69

   "global_p_024_002"  =>  [    2.0945,  1.0 ],  # GynCycle Param 70
   "global_p_024_003"  =>  [    9.28  ,  1.0 ],  # GynCycle Param 71
   "global_p_024_004"  =>  [ 3480.27  ,  1.0 ],  # GynCycle Param 72
   "global_p_024_005"  =>  [    0.972 ,  1.0 ],  # GynCycle Param 73
   "global_p_024_006"  =>  [ 1713.71  ,  1.0 ],  # GynCycle Param 74
   "global_p_024_007"  =>  [ 8675.14  ,  1.0 ],  # GynCycle Param 75
   "global_p_024_008"  =>  [    5.235 ,  1.0 ]   # GynCycle Param 76
}
pIniGuess = {
   "global_p_019_001"  =>  [    0.958 ,  1.0 ],  # GynCycle Param 61
   "global_p_020_001"  =>  [    0.925 ,  1.0 ],  # GynCycle Param 62
   "global_p_021_001"  =>  [    0.7567,  1.0 ],  # GynCycle Param 63
   "global_p_022_001"  =>  [    0.61  ,  1.0 ],  # GynCycle Param 64
   "global_p_023_001"  =>  [    0.543 ,  1.0 ],  # GynCycle Param 65

#  "global_p_024_001"  =>  [   51.558 ,  1.0 ]   # GynCycle Param 69

   "global_p_024_002"  =>  [    2.0945,  1.0 ]   # GynCycle Param 70
} if 0 < 1
pIniGuess = {
   "global_p_019_001"  =>  [    0.958 ,  1.0 ]   # GynCycle Param 61
} if 0 < 0
# pIniGuess = {}

nPar = 0
pIniGuess.keys.each do |key| 
  # next unless model.pId.index(key) 
  nPar += 1 if (model.pId.index(key) or model.y0Id.index(key)) 
end
# pidx   = [   2,   1 ]
# guess  = [ 3.5, 0.3 ]
# pscal  = [ 1.0, 1.0 ]
# nPar   = pidx.length

tspan = [] # (0..120).to_a  # [0.0,90.0]
0.upto(600) { |j| tspan << j/5.0 }
### tspan = [0.0,120.0]
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
if pidx==[] then
  pidx = (1..par.length).to_a
end

puts ""
puts "ModelDL: model.solve_var called with:"
puts "pidx : #{pidx}"
puts "par0 : #{model.par}"


tpoints, sol0 = model.solve_var tspan, y0, par, pidx

sout = File.open("rb_ModelDL_vareq_PAEON_V2_solution.dat","w")
model.save_current_solution sout
sout.close


spout = File.open("rb_ModelDL_vareq_PAEON_V2_soldp.dat","w")
pidx.each_with_index do |idx,j|
  psave = par[idx-1]
  par[idx-1] += $delp*psave
  tp, sol = model.solve_ode tspan, y0, par
  par[idx-1] = psave
  spout.printf("# dp = %.2f * %e = %e\n", $delp, psave, $delp*psave)
  # model.save_current_solution spout,idx
  spout.printf("# == %s ==\n", model.pId[idx-1])
  spout.printf("%-12s","Timepoint")
  model.yId.each do |label|
     label = "***" if label.length == 0
     spout.printf("\t%-12s", label)
  end
  spout.printf("\n")
  tp.each_with_index do |t,ell|
     if tp[ell] == tpoints[ell] then
       spout.printf("%.6e", t)
       y = sol[t] 
       yref = sol0[t]
       y.each_with_index { |val,k| spout.printf("\t%.6e", (y[k]-yref[k])/psave) }
       spout.printf("\n")
     end
  end
  spout.printf("\n\n")
end
spout.close


if $gsl_avail
   nspe = y0.length
   npar = pidx.length
   speout = File.open("rb_ModelDL_vareq_PAEON_V2_svd_species.dat","w")
   parout = File.open("rb_ModelDL_vareq_PAEON_V2_svd_parameter.dat","w")
   speout.printf("%-12s", "Timepoint")
   model.yId.each_with_index do |label,dmy|
      label = "***" if label.length == 0
      speout.printf("\tn%-11s\tx%-11s\ty%-11s\tz%-11s", 
                     label, label, label, label)
   end
   speout.printf("\n")
   parout.printf("%-12s", "Timepoint")
   pidx.each do |idx|
      label = "***"
      label = model.pId[ idx-1 ] if idx > 0
      label = "***" if label.length == 0
      parout.printf("\tn%-11s\tx%-11s\ty%-11s\tz%-11s", 
                     label, label, label, label)
   end
   parout.printf("\n")

   tpoints[1..-1].each do |t|

      mat = GSL::Matrix[ sol0[t][nspe..-1], nspe, npar ]

      if nspe < npar then
         v, u, s = mat.trans.SV_decomp
         # su = GSL::Matrix.diag(s)*u.trans
         # sv = v*GSL::Matrix.diag(s)
      else
         u, v, s = mat.SV_decomp
      end
      su = GSL::Matrix.diag(s)*u.trans
      sv = v*GSL::Matrix.diag(s)
      total = s.sum
      puts "  #{'%6.2f  %dx%d  %dx%d  %8.2e  %8.2e  %8.2e ... sum %8.2e' % 
              [t, u.shape, v.shape, s[0..2].to_a, total].flatten}"
      puts "  #{'%6.2f  %dx%d  %dx%d    %5.2f%%    %5.2f%%    %5.2f%% ...' % 
              [t, u.shape, v.shape, s[0]*100.0/total, (s[0]+s[1])*100.0/total, 
                                    (s[0]+s[1]+s[2])*100.0/total].flatten}"
      #
      speout.printf("% .6e", t) 
      0.upto(nspe-1) do |j| 
         speout.printf("\t% .6e\t% .6e\t% .6e\t% .6e",
                        su.col(j).norm, 
                        u[j,0]*s[0], 
                        u[j,1]*s[1], 
                        u[j,2]*s[2])  # (almost) principle components
      end
      speout.printf("\n")
      #
      parout.printf("% .6e", t)
      0.upto(npar-1) do |j|
         parout.printf("\t% .6e\t% .6e\t% .6e\t% .6e",
                        sv.row(j).norm, 
                        v[j,0]*s[0], 
                        v[j,1]*s[1], 
                        v[j,2]*s[2])
      end
      parout.printf("\n")
      #
   end
   speout.close
   parout.close
end

# ---------------------------------------------------------------------
# Result Output

#
#  already done; see above!
#
#sout = File.open("rb_ModelDL_vareq_PAEON_V2_solution.dat","w")
#model.save_current_solution sout
#sout.close

puts " "
puts "#{model.version}"
puts " "

