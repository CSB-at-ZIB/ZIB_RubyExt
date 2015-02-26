#! /usr/bin/env ruby

# =====================================================================

$LOAD_PATH << File.dirname(__FILE__)

require_relative '../../lib/ModelDL'
require_relative '../../lib/Experiment'
require_relative '../../lib/SysBioFit'

# =====================================================================

rtol = atol = 1.0e-4
rtol = atol = Float(ARGV[0]) if ARGV.length > 0
atol = Float(ARGV[1]) if ARGV.length > 1

# ---------------------------------------------------------------------
# Model/dynamic-load ODE

model = ModelDL.new 
model.t0 = -42.5
model.hmax = 0.0
model.inistep = 1.0e-6
model.rtol = rtol
model.atol = atol
# model.monitor = 1
# puts "#{model.inspect}"
puts "#{model.version}" # see at the end of this script


start = Time.now
t, sol = model.solve_ode [0.0, 150.0]
stop = Time.now

puts
puts "    #steps: #{t.length} (ODE solution in interval t=[0,150])"
puts "tolerances: #{model.tolerance}"
puts
puts "Time elapsed #{(stop - start)*1.0}s "
puts

