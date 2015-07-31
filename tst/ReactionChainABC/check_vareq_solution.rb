#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/Model'

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

def jacobian(t,y,par,pidx)

  dy = abc(t,y,par) # y is virtually too long, but this does not matter! 

  nspe = dy.length
  npar = par.length
  k1, k2 = par

  drdy = [ [  k1 , 0.0 , 0.0 ] ,
           [ 0.0 ,  k2 , 0.0 ] ]

  drdp = [ [ y[0] ,  0.0 ] ,
           [  0.0 , y[1] ] ]

  s = y[nspe..-1]  # m = 3 species (y[0], y[1], y[2]), and q = 2 parameters
                   #  ==>  s ~ y[m = 3], ..., y[(m+m*q-1) = (3 + 3*2 - 1) = 8]

  fy = [ [ -drdy[0][0]              , -drdy[0][1]              ,  drdy[0][2]              ] ,
         [  drdy[0][0] - drdy[1][0] ,  drdy[0][1] - drdy[1][1] ,  drdy[0][2] - drdy[1][2] ] ,
         [               drdy[1][0] ,               drdy[1][1] ,               drdy[1][2] ] ]

  fp = [ [ -drdp[0][0]              , -drdp[0][1]              ] ,
         [  drdp[0][0] - drdp[1][0] ,  drdp[0][1] - drdp[1][1] ] ,
         [               drdp[1][0] ,               drdp[1][1] ] ]

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

  [ dy ,  dS ].flatten!

end

# ---------------------------------------------------------------------
# Model/ODE setup

initvals = {
                   t0:  0.0 , 
                   y0: [ 1.0  , 0.0  , 0.0  ] ,
              y0label: [ "A0" , "B0" , "C0" ] ,

                  par: [  2.0 ,  1.0  ] ,
               plabel: [  "k1",  "k2" ] ,

                  jac: :jacobian
           }

model = Model.new :abc, initvals    # t0, y0, par, plabel
model.hmax = 0.0
model.inistep = 1.0e-4
model.rtol = $rtol
model.atol = $atol

# puts "#{model.version}"

# ---------------------------------------------------------------------
# Parameter Estimation/Identification

pIniGuess = {
               "k2" => [ 3.5 , 0.0 ],
               "k1" => [ 0.3 , 0.0 ]
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

sout = File.open("rb_vareq_solution.dat","w")
model.save_current_solution sout
sout.close

puts "#{model.version}"
puts " "

