require 'csv'

class Experiment

  attr_accessor :mtime, :mdata, :mweight
  attr_accessor :mlabel, :measmatrix
  attr_accessor :alpha

  def initialize( mtime = [], mdata = [], mweight = [], 
                  mlabel =[], measmatrix = nil, alpha = 0.95 )
    @mtime = mtime
    @mdata = mdata
    @mweight = mweight
    @mlabel = mlabel
    @measmatrix = measmatrix
    @alpha = alpha
  end


  def load_csv(csvfile)
    
    table = []

    if csvfile.is_a?(String) then
      table = CSV.read(csvfile, col_sep: ' ', converters: :numeric)
    elsif csvfile.is_a?(CSV) then
      csvfile.foreach do |row|
        table << row
      end
      # csvfile.close
    end

    table.each_with_index do |row,k|
      if k == 0 then
         row.shift
         @mlabel = row.map {|x| x.to_s }
      else
        @mtime << row.shift
        @mdata << row
        weights = []
        row.each { |dummy| weights << 1.0 }
        @mweight << weights
      end
    end

    # @mdata.flatten!
    #
    # @mweight = []
    # @mdata.each { |dummy| @mweight << 1.0 }
  end


  def load_data(csvfile)
    
    table = []

    if csvfile.is_a?(String) then
      table = CSV.read(csvfile, col_sep: ' ', converters: :numeric)
    elsif csvfile.is_a?(CSV) then
      csvfile.foreach do |row|
        table << row
      end
      # csvfile.close
    end

    @mlabel = []
    @mtime = []
    @mdata = []
    @mweight = []

    table.each_with_index do |row,k|
       if k == 0 then
          row.shift
          while row.length > 1 do
             @mlabel << row.shift
             row.shift
          end
       else
          @mtime << row.shift if row.length > 0

          measurements = []
          weights = []
          while row.length > 1 do
             measurements << row.shift
             weights << row.shift
          end
          @mdata << measurements unless measurements==[]
          @mweight << weights unless weights==[]
       end
    end

    # @mdata.flatten!
    # @mweight.flatten!
  end

end
