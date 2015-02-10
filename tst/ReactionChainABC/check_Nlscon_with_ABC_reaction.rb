#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/Model'
require_relative '../../lib/Experiment'
require_relative '../../lib/SysBioFit'

# ---------------------------------------------------------------------
# ODE system

def abc(t,y,par)

  k1,k2 = par

  r1 = k1 * y[0]
  r2 = k2 * y[1]

  # dy0     dy1      dy2
  [ -r1 , +r1 - r2 , +r2 ]

end

#

def jacobian(t,y,par)

  q = par.length
  k1, k2 = par


  drdy = [ [  k1 , 0.0 , 0.0 ] ,
           [ 0.0 ,  k2 , 0.0 ] ]

  drdp = [ [ y[0] ,  0.0 ] ,
           [  0.0 , y[1] ] ]

  s = y[3..-1]  # 3 species (y[0], y[1], y[2]), and 2 parameters
                #  ==>  s ~ y[3], ... , y[8]

  fy = [ [ -drdy[0][0]              , -drdy[0][1]              ,  drdy[0][2]              ] ,
         [  drdy[0][0] - drdy[1][0] ,  drdy[0][1] - drdy[1][1] ,  drdy[0][2] - drdy[1][2] ] ,
         [               drdy[1][0] ,               drdy[1][1] ,               drdy[1][2] ] ]

  fp = [ [ -drdp[0][0]              , -drdp[0][1]              ] ,
         [  drdp[0][0] - drdp[1][0] ,  drdp[0][1] - drdp[1][1] ] ,
         [               drdp[1][0] ,               drdp[1][1] ] ]

  dy = abc(t,y,par)

  dS = []
  q.times do |ell|
    3.times do |j|
      sum = fp[j][ell]
      3.times {|k| sum += fy[j][k]*s[q*ell + k] }
      dS[q*ell + j] = sum
    end
  end

  [ dy ,  dS ].flatten!

end

# ---------------------------------------------------------------------
# Model/ODE setup

initvals = {
                   t0:  0.0 , 
                   y0: [ 1.0  , 0.0  , 0.0  ] ,
              y0label: [ "A0" , "B0" , "C0" ],

                  par: [  2.0 ,  1.0  ] ,
               plabel: [  "k1",  "k2" ]
           }

model = Model.new :abc, initvals    # t0, y0, par, plabel

# puts "#{model.version}"

# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_ABC_data.dat"

data = Experiment.new 
data.load_csv fname

# ---------------------------------------------------------------------
# Parameter Estimation/Identification

pIniGuess = {
               "k2" => [ 3.5 , 1.0 ],
               "k1" => [ 0.3 , 1.0 ]
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
fit.rtol = 1.0e-5
fit.printlevel = 2
fit.pfname = "rb_Nlscon_with_ABC_parameter.dat"
fit.sfname = "rb_Nlscon_with_ABC_solution.dat"
# fit.nonlin = 2
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

