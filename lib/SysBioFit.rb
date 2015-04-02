$LOAD_PATH << File.dirname(__FILE__)

require_relative '../ext/Nlscon/Nlscon'
# require_relative './Model'
# require_relative './Experiment'

# ---------------------------------------------------------------------
# ---------------------------------------------------------------------

class SysBioFit < Nlscon

  attr_accessor :pidx, :pguess, :pscal, :parub, :parlb
  attr_accessor :nitmax, :nonlin, :arrhenius, :printlevel, :jacgen
  attr_accessor :nitmodulo, :sfname, :pfname

  def initialize( ary, pidx = [], pguess = [], pscal = [], 
                  parub = [], parlb = [], ptol = 1.0e-3, 
                  nitmax = 35, nonlin = 3, jacgen = 3, printlevel = 0)
    super(ary)
    @pidx = pidx
    @pguess = pguess
    @pscal = pscal
    @parub = parub
    @parlb = parlb
    @nitmax = nitmax
    @nonlin = nonlin
    @jacgen = jacgen
    # @arrhenius = arrhenius
    @printlevel = printlevel
    @sfname = nil
    @pfname = nil
    @plabel = nil
    @model = nil
    @tspan = nil
    @mtime = nil
    @midx = nil
    @mdata = nil
    @nitmodulo = 5
    self.rtol = ptol
    self.iopt = { "mprmon" => 0, "mprerr" => 0, "mprsta" => 0, "qstat" => 1,
                  "qrank1" => 0, "nonlin" => 0, "jacgen" => 3, "qnscal" => 0 }
    self.iwk = { "nitmax" => 35 }
    # self.rwk = { "cond" => 1.0e+16 }
  end


  def save_parameter(fout,par,iter)
    return if fout.closed?

    if iter == 0 then
      fout.printf("%-12s", "Iteration")
      if @plabel != nil and 
         @plabel.compact.length == par.compact.length 
      then
        @plabel.each_with_index do |label,idx|
          label="***" if label.length == 0
          fout.printf("\t%-12s", label)
          fout.printf("%c", "+") if @pidx.count(idx+1) != 0
        end
      else
        par.each_with_index do |val,idx|
          fout.printf("\t%12d", idx+1)
          fout.printf("%c", "+") if @pidx.count(idx+1) != 0
        end
      end
      fout.printf("\t%-12s", "SKAP")
      fout.printf("\t%-12s", "DLEVF")
      fout.printf("\n")
    end
    
    fout.printf(" %6d     ", iter)
    par.each { |val| fout.printf("\t%.6e", val) }
    fout.printf("\t%.6e", (iter > 0) ? self.rwk["skap"] : 0.0)
    fout.printf("\t%.6e", (iter > 0) ? self.rwk["dlevf"] : 0.0)
    fout.printf("\n") 
  end


  def set_current_par(x)
    # idxMax = x.length
    y0 = @model.y0ode.clone unless @model==nil
     p = @model.par0.clone unless @model==nil
    @pidx.each_with_index do |j,idx| 
         p[j-1] = x[idx] if j > 0
       y0[-j-1] = x[idx] if j < 0
    end
    [y0, p]
  end


  def gn_fcn(n, m, mcon, x)

    # t = par = sol = y = nil
    #n.times {|j| @par[ @pidx[j]-1 ] = x[j] }
    y0, par = set_current_par(x)
    
    t, sol = @model.solve_ode(@tspan, y0, par)

    fsim = []

    return if sol==nil

    # (m-mcon).times do |k|
    @mtime.length.times do |k|
      next if @midx[k]==[]
      y = sol[ @mtime[k] ]
      @midx[k].each_with_index do |j,idx|
        # puts "sim = #{y[j]}       meas = #{@mdata[k][idx]}"
        fsim << y[j] unless j==nil
      end
    end

    fsim
    # fsim.flatten!
  end


  def set_task(task={})

    return unless task.has_key?(:model) or
                  task.has_key?(:data) or
                  task.has_key?(:guess)

    @model = task[:model]
    # @plabel = task[:guess].keys 
    @plabel = task[:model].pId
    @mtime = task[:data].mtime.map { |x| x.to_f }
    @mdata = task[:data].mdata
    @midx = []
    task[:data].mdata.each do |ary|
      row = []
      ary.each_with_index do |idx,j|
        next if idx.is_a?(String)
        label = task[:data].mlabel[j]
        k = task[:model].yId.index(label)
        row << k
      end
      @midx << row
    end
    @tspan = [task[:model].t0].concat(@mtime.sort).uniq

    if task[:guess].is_a?(Array) and task.has_key?(:pidx) then
      @pidx = task[:pidx]
      @pguess = task[:guess]
      @pscal = []
      if task.has_key?(:pscal) then
        @pscal = task[:pscal]
      end
    elsif task[:guess].is_a?(Hash) then
      @pidx = []
      @pguess = []
      @pscal = []
      task[:guess].each do |key,ary|
         j = task[:model].pId.index(key)
         k = task[:model].y0Id.index(key)
         if j then
           @pidx << j+1    # positive indices > 0 for parameters
         elsif k then
           @pidx << -k-1   # negative indices < 0 for initial conditions
         else
           next
         end
         @pguess << ary[0]
         @pscal << ary[1]
      end
    else
      return
    end

    nPar = @pidx.length
    # mTotal = @mdata.flatten.length - @mdata.flatten.count("n/a")
    # mFit = mTotal - task[:data].mweight.count(0.0)

    self.x = (@pguess.compact.length == nPar) ? @pguess : [].fill(0.0,0,nPar)
    self.xscal = (@pscal.length == nPar) ? @pscal : [].fill(0.0,0,nPar)

    measurements = task[:data].mdata.flatten
    measurements.delete_if { |val| val.is_a?(String) }
    weights = task[:data].mweight.flatten
    weights.delete_if { |val| val.is_a?(String) }

    self.fobs = measurements
    self.fscal = weights if weights.length == measurements.length

  end



  def identify_par(task={})
  
    return unless task.has_key?(:model) or 
                  task.has_key?(:data) or 
                  task.has_key?(:guess) 
    pout = nil
    sout = nil
    pout = File.open(@pfname, "w") if @pfname
    sout = File.open(@sfname, "w") if @sfname

    set_task(task)

    self.f = method(:gn_fcn)
    # self.df = method(:gn_jac)

    self.iopt = { "mprmon" => @printlevel, "mprstat" => 1, "mprerr" => [@printlevel,3].min, 
                  "jacgen" => @jacgen, "nonlin" => @nonlin, "qstat" => 1 }
    self.iwk = { "nitmax" => @nitmax }

    tinterval = [@tspan[0],@tspan[-1]]
    nitmax = self.iwk["nitmax"]
    iter = 0
    status = -1

    y0, par = set_current_par(self.x) if (pout or sout)
    save_parameter(pout,par,iter) if pout
    @model.solve_ode(tinterval,y0,par) if sout
    @model.save_current_solution(sout,iter) if sout
    
    while status == -1 && iter < nitmax do
      iter += 1
      status = self.iterate   # { |x| gn_fcn(x) }

      # puts "identify_par: status = #{status}"

      y0, par = set_current_par(self.x) if (sout or pout)
      save_parameter(pout, par, iter) if pout
      if sout and iter.modulo(@nitmodulo) == 0.0 then
        @model.solve_ode( tinterval, y0, par)
        @model.save_current_solution(sout, iter)
      end

      printf("%5d\n",iter) if @printlevel > 0
    end

    y0, par = set_current_par(self.x)
    @model.solve_ode( tinterval, y0, par )
    @model.save_current_solution(sout, -1) if sout

    pout.close if pout
    sout.close if sout

    { "par" => self.x , "iwk" => self.iwk , "rwk" => self.rwk , "iopt" => self.iopt ,
      "pidx" => @pidx, "t" => @model.solution.keys , "y" => @model.solution.values }
  end


  def dsign(x,y)
    return   x.abs  if y.to_f >= 0.0
    return -(x.abs) if y.to_f < 0.0
  end


  def ncjac(y)

    ajmin = self.rwk["ajmin"]
    ajdel = self.rwk["ajdel"]

    n = self.y.length 

    z0 = []
    z0 = gn_fcn(n,0,0,y)

    return if z.compact == nil

    aa = []
    y.each_with_index do |w,k|

       u = dsign( ajdel*([w.abs,ajmin,self.xscal[k]].max), w )
       y[k] = w + u
     
       z = []
       z = gn_fcn(n,0,0,y)

       return if z.compact == nil or z.length != z0.length
 
       y[k] = w
 
       vec = []
       z.each_with_index do |v,j|
         vec << (v - z0[j]).to_f / u
       end

       aa << vec

    end

    aa

  end


  def ncjcf(y)

    etaini = ( self.rwk["etaini"] > 0.0 ) ? self.rwk["etaini"] : 1.0e-6
    etadif = ( self.rwk["etadif"] > 0.0 ) ? self.rwk["etadif"] : 1.0e-6
    etamax = Math.sqrt( Math.sqrt(10.0*Float::EPSILON) )
    etamin = etamax * Math.sqrt( 10.0*Float::EPSILON )
    eta = self.rwk["eta"] || [].fill(etaini,0,y.length)

    n = y.length

    z0 = []
    z0 = gn_fcn(n,0,0,y)

    aa = []
    y.each_with_index do |w,k|

       is = 0
       qfine = false

       while not qfine do

          w = y[k]
          u = dsign( eta[k]*self.xscal[k], w )

          y[k] = w + u

          z = []
          z = gn_fcn(n,0,0,y)

          return if z.compact == nil or z.length != z0.length

          y[k] = w

          sum = 0.0
          vec = []
          z.each_with_index do |v,j|
             hg  = [ v.abs, z0[j].abs ].max
            fhj  = v - z0[j]
            sum += (fhj/hg)**2 if hg != 0.0
            vec << fhj / u
          end

          sum = Math.sqrt( sum / z.length )
          qfine = true
          if sum != 0.0 and is == 0 then
            eta[k] = [ etamax, [etamin, Math.sqrt(etadif/sum)*eta[k]].max ].min
            is = 1
            # qfine = (conv < small2) or (sum >= etamin)
            qfine = ( sum >= etamin )
          end

       end

       aa << vec

    end  #  y.each_with_index

    aa

  end

  
  def compute_sensitivity(task={}, nrm="l2")
  
    return unless task.has_key?(:model) or 
                  task.has_key?(:data) or 
                  task.has_key?(:guess) 
    # tout = nil
    # tout = File.open(@sensfname, "w") if @sensfname

    set_task(task)

    aa = []
    jcnrm = []

    if self.iopt["jacgen"] == 3 then
      aa = ncjcf(self.x.clone)
    elsif self.iopt["jacgen"] == 2 then
      aa = ncjac(self.x.clone)
    end

    if nrm == "linf" then
      aa.each do |ary|
        linf = ary.map {|val| val.abs }.max
        ary.map! {|val| val/linf } if linf > 0.0
        jcnrm << linf
      end
    elsif nrm == "l1" then
      aa.each do |ary|
        l1  = ary.map {|val| val.abs }.inject(:+)
        l1 /= ary.length unless ary.length == 0
        ary.map! {|val| val/l1 } if l1 > 0.0
        jcnrm << l1
      end
    else
      nrm = "l2"
      aa.each do |ary|
        l2 = Math.sqrt( ary.map {|val| (val.abs)**2 }.inject(:+) )
        l2 /= Math.sqrt( ary.length ) unless ary.length == 0
        ary.map! {|val| val/l2 } if l2 > 0.0
        jcnrm << l2
      end
    end

    { "pidx" => @pidx, "jacobian" => aa, "jcnrm" => jcnrm, "nrm" => nrm }

  end

end

