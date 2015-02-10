#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/Model'
require_relative '../../lib/Experiment'
require_relative '../../lib/SysBioFit'

# ---------------------------------------------------------------------
# ODE system

def predator_prey(t,y,par)

  alph,bet,gam,del = par

  n1 = y[0]
  n2 = y[1]

  #     d n1 / dt             d n2 / dt
  [ n1*(alph - bet*n2) , -n2*(gam - del*n1) ]

end

# ---------------------------------------------------------------------
# Model/ODE setup

initvals = {
                    t0:     1900.0 ,  

                    y0: [  30.0 ,   4.0  ] ,
               y0label: [ "n1_0", "n2_0" ] ,

                   par: [    0.5 ,  0.02 ,    1.0 ,   0.02  ] ,
                plabel: [ "alpha", "beta", "gamma", "delta" ]
           }

model = Model.new :predator_prey, initvals   # t0, y0, par, plabel

# puts "#{model.inspect}"
# exit


# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "hare_lynx_data.txt"

data = Experiment.new 
data.load_csv fname

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

# puts "#{data.inspect}"

# ---------------------------------------------------------------------
# Parameter Estimation/Identification

pIniGuess = {
               "alpha" => [ 0.5  ,  1.0 ],  #  0.1  ],
                "beta" => [ 0.02 ,  1.0 ],  #  0.01 ],
               "gamma" => [ 1.0  ,  1.0 ],  #  1.0  ],
               "delta" => [ 0.02 ,  1.0 ]  #  0.01 ] 
             #   "n1_0" => [ 30.0 ,  1.0 ],  # 10.0  ] 
             #   "n2_0" => [  4.0 ,  1.0 ]   #  1.0  ]
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

nlscon = SysBioFit.new [nPar,mTotal,mFit]
nlscon.rtol = 1.0e-4
nlscon.printlevel = 5 # 3 # 2
nlscon.pfname = "rb_Nlscon_with_PredPrey_parameter.dat"
nlscon.sfname = "rb_Nlscon_with_PredPrey_solution.dat"
# nlscon.nonlin = 4
# nlscon.rwk = { "cond" => 1.0e+7 } 

# puts "nlscon.rwk = #{nlscon.rwk}"

id_task = { model: model, data: data,  guess: pIniGuess }
           # pidx: pidx, guess: guess, pscal: pscal }
         
estim = nlscon.identify_par id_task 

# ---------------------------------------------------------------------
# Result Output

printf "\n"
printf "   %12s     %12s   %12s     %12s\n", 
       "Parameter", "confid_lower", "*p_estim*", "confid_upper"
printf "   %12s     %12s   %12s     %12s\n", 
       "---------", "------------", "---------", "------------"

estim["pidx"].each_with_index do |j,idx|
  label = (j > 0) ? model.pId[j-1] : model.y0Id[-j-1] 
  printf "  %12s     % 12.6f    % 12.6f    % 12.6f\n", 
          label[0..11], 
          estim["rwk"]["xl"][idx], 
          estim["par"][idx], 
          estim["rwk"]["xr"][idx]
end

printf "\n"

# puts " "
# puts "                   n1_0         n2_0"
# puts "                ----------   ----------"
# puts " confid_lower:  #{'% .6f    % .6f' % estim["rwk"]["xl"][-2..-1]}"
# puts "  *p_estim*  :  #{'% .6f    % .6f' % estim["par"][-2..-1]}"
# puts " confid_upper:  #{'% .6f    % .6f' % estim["rwk"]["xr"][-2..-1]}"
# puts " "
# puts "                  alpha         beta        gamma        delta"
# puts "                ----------   ----------   ----------   ----------"
# puts " confid_lower:  #{'% .6f    % .6f    % .6f    % .6f' % estim["rwk"]["xl"]}"
# puts "  *p_estim*  :  #{'% .6f    % .6f    % .6f    % .6f' % estim["par"]}"
# puts " confid_upper:  #{'% .6f    % .6f    % .6f    % .6f' % estim["rwk"]["xr"]}"
# puts " "
puts "incomp. kappa:    #{'% .4e' % estim["rwk"]["skap"]}"
puts "achieved rtol:    #{'% .4e' % estim["rwk"]["prec"]}"
puts " "

# puts "#{estim["iopt"]}"
# puts "#{estim["rwk"]}"

# puts "#{estim["t"]}"
# puts "#{estim["y"]}"

# ---------------------------------------------------------------------
# Sensitivity Matrix Output (w.r.t. initial guess)

sens = nlscon.compute_l2sensitivity id_task, "l2"


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
end

channel.shift
channel.each { |out| out.close }

