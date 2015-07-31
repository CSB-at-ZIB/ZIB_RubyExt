$LOAD_PATH << File.dirname(__FILE__)

require_relative '../ext/LimexDL/LimexDL'

# ---------------------------------------------------------------------

class ModelDL < LimexDL

  attr_accessor :t0, :y0ode, :y0Id, :par0
  # attr_accessor :par, :rtol, :atol, :inistep, :hmax
  # attr_accessor :y0, :rtol, :atol, :inistep, :hmax

  def initialize( initvals = {},
                  rtol = 1.0e-9, atol = 1.0e-9, 
                  inistep = 0.0, maxstep = 0.0 )

    super()
    y0 = initvals.has_key?(:y0) ? initvals[:y0] : []
    par = initvals.has_key?(:par) ? initvals[:par] : []
    t0 = initvals.has_key?(:t0) ? initvals[:t0] : 0.0

    @pidx = []
    @y0ode = self.y0
    @y0Id = self.yId
    if initvals.has_key?(:y0label) and
       y0.length == initvals[:y0label].length
    then
      @y0Id = initvals[:y0label]
    end
    @t0 = t0
    @par0 = self.par
    # @rtol = rtol
    # @atol = atol
    # @inistep = inistep
    # @maxstep = maxstep
    # @y0 = self.y0
    self.par = par if par.compact != []
    self.y0 = y0 if y0.compact != []
    self.rtol = rtol
    self.atol = atol
    self.inistep = inistep
    self.hmax = maxstep
  end


  def solve_var(tspan, y0=[], par=[], pidx=[])
  #
    @pidx = []
    self.y0 = y0 if y0.compact != []
    self.par = par if par.compact != []
    return [nil,nil] if pidx.compact == []
    ifail = self.srun( tspan.sort, pidx ) 
    return [nil,nil] unless ifail.is_a?(Array) and ifail[0] == 0
    @pidx = pidx
    [ self.steps, self.solution ]
  #
  end


  def solve_ode(tspan, y0=[], par=[])
  #
    @pidx = []
    self.y0 = y0 if y0.compact != []
    self.par = par if par.compact != []
    ifail = self.run( tspan.sort ) 
    return [nil,nil] unless ifail.is_a?(Array) and ifail[0] == 0
    [ self.steps, self.solution ]
  #
  end


  def save_current_solution(fout, iter=-2)
  #
    return if fout.closed?

    fout.printf("# == Iter %4d ==\n", iter) if iter != -2
    fout.printf("%-12s","Timepoint")
    self.yId.each do |label|
       label = "***" if label.length == 0
       fout.printf("\t%-12s", label)
    end
    if @pidx.compact != [] and 
       self.solution.first[1].length == self.yId.length*(1+@pidx.length) 
    then
       @pidx.each_with_index do |ell,dmy|
          pId = "***"
          pId = self.pId[ ell-1 ] if ell > 0
          self.yId.each do |label|
             label = "***" if label.length == 0
             str = "#{'%5s/%5s' % [label, pId]}"
             fout.printf("\t%-12s", str)
          end
       end
    end
    fout.printf("\n")
    # fout.printf("%.6e", self.interval[0])
    # self.y0.each { |val| fout.printf("\t%.6e", val) }
    # fout.printf("\n")
    self.steps.each do |t|
       fout.printf("%.6e", t)
       y = self.solution[t]
       y.each { |val| fout.printf("\t%.6e", val) }
       fout.printf("\n")
    end
    fout.printf("\n\n")
  #
  end

end
