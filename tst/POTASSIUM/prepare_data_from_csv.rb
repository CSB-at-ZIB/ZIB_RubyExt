#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require 'csv'
require_relative '../../lib/ModelDL'

# ---------------------------------------------------------------------

$dpath = "."

$fname = "#{$dpath}/kuh_a_all.csv"
$data = { 
             "Timepoint"  =>  "Time" ,
             "K_ECF"      =>  "K_ECF" ,
             "K_ICF"      =>  "K_ICF",
             "s29"        =>  "Glucose" ,
             "Insulin"    =>  "Insulin",
             "K_urin"     =>  "K_urin"
        }

$tspan = [5.0, 70.0]

$dout = "rb_prepare_data_data.dat"

# ---------------------------------------------------------------------

$d = {}

  csv = CSV.open($fname, col_sep: "\t", headers: true)
  csv.to_a.map do |row|
     h = row.to_hash 
     puts "#{h}"
     t = Float(h[$data["Timepoint"]]) rescue nil
     next unless t
     if not $d.has_key?(t) then
       $d[t] = {}
     end
     $data.each_pair do |spe,label|
        if spe != "Timepoint" then
           y = Float(h[label]) rescue nil
           if $d[t].has_key?(spe) then
             $d[t][spe] << y
           else
             $d[t][spe] = [y]
           end
        end
     end
  end

  csv.close


$d.each_key do |t|
  m = 0
  $d[t].each_value { |ary| if ary.length > m then m = ary.length end }
  $d[t]["_arylength_"] = m
end

puts "#{$d}"

# ---------------------------------------------------------------------
# Write the csv data (in a format supported by the class "Experiment")

File.open($dout, "w") do |fout|
  fout.printf("%-12s", "Timepoint")
  $data.each_key do |spe| 
    fout.printf("\t%-12s\t%-12s", spe, "SD") unless spe == "Timepoint" 
  end
  fout.printf("\n")
  $d.each_key do |t|
    tp = Float(t) rescue nil
    next unless tp
    1.upto($d[t]["_arylength_"]) do |k|
      fout.printf("%.6e", tp)
      $data.each_key do |spe|
        next if spe == "Timepoint"
        ary = $d[t][spe]
        val = ary[k-1]
        if val then
          fout.printf("\t%.6e\t%.6e", val, 1.0)
        else
          fout.printf("\t%-12s\t%-12s", "n/a", "n/a")
        end
      end
      fout.printf("\n")
    end
  end
  fout.printf("\n")
end

# ---------------------------------------------------------------------
# Load the dynamic load model (here, in this case, a POTASSIUM model!)

model = ModelDL.new 

model.solve_ode($tspan)

fout = File.open("rb_prepare_data_solution.dat", "w")
model.save_current_solution(fout)
fout.close

