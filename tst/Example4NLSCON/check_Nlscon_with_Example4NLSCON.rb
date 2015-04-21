#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

# require_relative '../../lib/Model'
require_relative '../../lib/Experiment'
# require_relative '../../lib/SysBioFit'
require_relative '../../ext/Nlscon/Nlscon'

# ---------------------------------------------------------------------
# FCN

def f(n,m,mcon,x)
  fx = []
  x1, x2, x3 = x

  fx[0] = x1 - x2**2 - x3**2 + 100.0

  m.times do |j|
    h = (2.0*( (j+1) - 1.0 ) - x3) / x2
    fx[j+1] = x1*Math::exp(-h*h/2.0)
  end

  fx
end

#

$rtol = 1.0e-8

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
$eta = Math::sqrt(10.0*$rtol) if $eta <= 0.0

# ---------------------------------------------------------------------
# Measurement/Experiment Data

fname = "rb_Nlscon_with_Example4NLSCON_data.dat"

data = Experiment.new 
### data.load_csv  fname  # without SD
data.load_data fname  # with SD

$measurements = data.mdata.flatten
# $fweights = [].fill(0.0,0,$measurements.length)
$fweights = data.mweight.flatten

# ---------------------------------------------------------------------
# Parameter Estimation/Identification

$x     = [ 1.0 , 2.0 , 5.0 ]
$xscal = [ 0.0 , 0.0 , 0.0 ]

$pfname = "rb_Nlscon_with_Example4NLSCON_parameter.dat"
$sfname = "rb_Nlscon_with_Example4NLSCON_solution.dat"

$nPar   = 3
$mTotal = 1 + $measurements.length
$mFit   = $measurements.length

fit = Nlscon.new [$nPar,$mTotal,$mFit]
fit.f = method(:f)
fit.rtol = $ptol
fit.iopt = {
              "mprerr" => 3, "mprmon" => 3, "mprstat" => 1,
              "jacgen" => 2, "nonlin" => 3, "qstat" => 1,
              "qnscal" => 0, "qrank1" => 0 
           }
fit.iwk = { "nitmax" => 195 }
fit.rwk = { "ajdel" => $eta }
# fit.rwk = { "cond" => 1.0e+16, "ajdel" => $eta }

fit.x = $x
fit.xscal = $xscal
fit.fobs = $measurements
fit.fscal = $fweights

nitmx = fit.iwk["nitmax"]
status = -1
iter = 0

while status == -1 && iter < nitmx do

  iter += 1
  status = fit.iterate

end

# ---------------------------------------------------------------------
# Result Output

