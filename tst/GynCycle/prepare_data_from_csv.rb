#! /usr/bin/env ruby

$LOAD_PATH << File.dirname(__FILE__)

require 'csv'
require_relative '../../lib/ModelDL'

# ---------------------------------------------------------------------

$dpath = "./Data/PharmaCompany"

$data = { 
            "E2"         =>  "#{$dpath}/Hormones-Estradiol.csv" ,
            "P4"         =>  "#{$dpath}/Hormones-Progesterone.csv" ,
            "LH_blood"   =>  "#{$dpath}/Hormones-LH.csv" ,
            "FSH_blood"  =>  "#{$dpath}/Hormones-FSH.csv"
        }

$tspan = [-42.5, 50.0]

$dout = "rb_Nlscon_with_PAEON_V2_data.dat"

# ---------------------------------------------------------------------

$d = {}

$data.each_pair do |spe,fname|

  csv = CSV.open(fname, col_sep: ',', headers: true)
  csv.to_a.map do |row| 
     x = false
     x = true if Integer(row[0]) rescue false
     break unless x
     row.to_hash.each_pair do |t,val|
       next unless t
       if spe == "E2" then     # silly unit change in Brigitte's data for E2
         y = 272.39*Float(val)/1000.0 rescue nil
       elsif spe == "P4" then  # ... and P4 <sigh>
         y = 314.47*Float(val)/1000.0 rescue nil
       else
         y = Float(val) rescue nil
       end
       # next unless y
       if $d.has_key?(t) then
         if $d[t].has_key?(spe) then
           $d[t][spe] << y
         else
           $d[t][spe] = [y] 
         end
       else
         $d[t] = {}
	 $d[t][spe] = [y]
       end
     end
  end
  csv.close

end 

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
          fout.printf("\t%.6e\t%.6e", val, val)  # , 1.0)
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
# Load the dynamic load model (here, in this case, a PAEON_V2 model!)

model = ModelDL.new

model.solve_ode($tspan)

fout = File.open("rb_prepare_data_solution.dat", "w")
model.save_current_solution(fout)
fout.close

puts "#{model.version}"

