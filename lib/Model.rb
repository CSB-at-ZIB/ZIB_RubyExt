$LOAD_PATH << File.dirname(__FILE__)

require_relative '../ext/Limex/Limex'

# ---------------------------------------------------------------------

class Model < Limex

  attr_accessor :t0, :y0ode, :y0Id, :par0
  attr_accessor :yId, :pId, :par, :version
  # attr_accessor :rtol, :atol, :inistep, :hmax
  # attr_accessor :y0, :rtol, :atol, :inistep, :hmax

  def initialize( odefcn = nil, initvals = {},
                  rtol = 1.0e-9, atol = 1.0e-9, 
                  inistep = 0.0, maxstep = 0.0 )

    now = Time.now()
    modelname = odefcn.is_a?(Symbol) ? odefcn.to_s : "n/a" 
    y0 = initvals.has_key?(:y0) ? initvals[:y0] : []
    par = initvals.has_key?(:par) ? initvals[:par] : []

    super(y0)
    @pidx = []
    @odejac = nil
    @odefcn = odefcn.is_a?(Symbol) ? method(odefcn) : nil
    @odefcn = odefcn if odefcn.is_a?(Proc)
    @y0ode = y0
    @version = ["#{modelname}", "#{now.asctime}", "#{now.to_i}", "n/a"]
    if initvals.has_key?(:version) and
        initvals[:version].length == 4
    then
      @version = initvals[:version]
    end

    @yId = []
    y0.each_with_index { |x,j| @yId << (j+1).to_s }

    if initvals.has_key?(:y0label) and
        y0.length == initvals[:y0label].length
    then
      @y0Id = initvals[:y0label]
    else
      @y0Id = []
      y0.each_with_index { |x,j| @y0Id << (j+1).to_s }
    end

    @par0 = @par = par
    if initvals.has_key?(:plabel) and
       par.length == initvals[:plabel].length 
    then
      @pId = initvals[:plabel]
    else
      @pId = []
      par.each_with_index { |x,j| @pId << (j+1).to_s }
    end

    @t0 = initvals.has_key?(:t0) ? initvals[:t0] : 0.0

    if initvals.has_key?(:jac) then
      @odejac = initvals[:jac].is_a?(Symbol) ? method(initvals[:jac]) : nil
      @odejac = initvals[:jac] if initvals[:jac].is_a?(Proc)
    end

    if @odejac.nil? then
      @version << "no vareq"
    else
      @version << "with vareq"
    end

    # @rtol = rtol
    # @atol = atol
    # @inistep = inistep
    # @maxstep = maxstep
    # self.y0 = y0
    self.rtol = rtol
    self.atol = atol
    self.inistep = inistep
    self.hmax = maxstep
  end


  def odefcn=(fcn)
    @odefcn = fcn.is_a?(Symbol) ? method(fcn) : nil
    @odefcn = fcn if fcn.is_a?(Proc)
  end

  def odejac=(jac)
    @odejac = jac.is_a?(Symbol) ? method(jac) : nil
    @odejac = jac if jac.is_a?(Proc)
  end


  def solve_var(tspan, y0=[], par=[], pidx=[])
  #
    @pidx = []
    self.y0 = y0 if y0.compact != []
    self.par = par if par.compact != []
    return [nil,nil] if @odejac.nil? or pidx.compact == []
    ifail = self.srun( tspan.sort, pidx ) do |t,y|
       if y.length == self.y0.length then
          @odefcn.call(t,y,self.par)
       else
          @odejac.call(t,y,self.par,pidx)
       end
    end
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
    return [nil,nil] if @odefcn.nil?
    ifail = self.run( tspan.sort ) do |t,y| 
        @odefcn.call(t,y,self.par)
    end
    return [nil,nil] unless ifail.is_a?(Array) and ifail[0] == 0
    [ self.steps, self.solution ]
  #
  end


  def save_current_solution(fout, iter=-2)
  #
    return if fout.closed?

    fout.printf("# == Iter %4d ==\n", iter) if iter != -2
    fout.printf("%-12s","Timepoint")
    self.y0.each_with_index do |val,idx| 
       fout.printf("\t%12d", idx+1)
    end
    if @pidx.compact != [] and 
       self.solution.first[1].length == self.y0.length*(1+@pidx.length) 
    then
       @pidx.each_with_index do |ell,dmy|
          self.y0.each_with_index do |val,idx|
             str = "#{'%05d/%05d' % [idx+1, ell]}"
             fout.printf("\t%12s", str)
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
