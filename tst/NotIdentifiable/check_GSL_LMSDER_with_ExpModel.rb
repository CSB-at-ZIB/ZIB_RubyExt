#! /usr/bin/env ruby

require("gsl")
include GSL::MultiFit

# x: Vector, list of the parameters to determine
# t, y, sigma: Vectors, observational data
# f: Vector, function to minimize
procf = Proc.new { |x, t, y, sigma, f|
  a = x[0]
  lamb = x[1]
  b = x[2]
  n = t.size
  for i in 0...n do
    yi = a*Math::exp(-lamb*t[i]) + b
    f[i] = (yi - y[i])/sigma[i]
  end
}

# jac: Matrix, Jacobian
procdf = Proc.new { |x, t, y, sigma, jac|
  a = x[0]
  lamb = x[1]
  n = t.size
  for i in 0...n do
    ti = t[i]
    si = sigma[i]
    ei = Math::exp(-lamb*ti)
    jac.set(i, 0, ei/si)
    jac.set(i, 1, -ti*a*ei/si)
    jac.set(i, 2, 1.0/si)
  end
}

f = Function_fdf.alloc(procf, procdf, 2)

n = 40

# Create data
r = GSL::Rng.alloc()
t = GSL::Vector.alloc(n)
y = GSL::Vector.alloc(n)
sigma = GSL::Vector.alloc(n)
# true parameters: { A: 5, lamb: 0.1, b = 1.0 }
for i in 0...n do
  t[i] = i
  y[i] = 1.0 + 5*Math::exp(-0.1*t[i]) + r.gaussian(0.1)
  sigma[i] = 0.1
end

f.set_data(t, y, sigma)
x = GSL::Vector.alloc(1.0, 0.0, 0.0)    # initial guess
np = 3

solver = FdfSolver.alloc(FdfSolver::LMSDER, n, np)

solver.set(f, x)

iter = 0
solver.print_state(iter)
begin
  iter += 1
  status = solver.iterate
  solver.print_state(iter)
  status = solver.test_delta(1e-4, 1e-4)
end while status == GSL::CONTINUE and iter < 500

covar = solver.covar(0.0)
position = solver.position
chi = solver.f.dnrm2
dof = n - np
printf("chi2/dof = %.6f\n", chi**2/dof)
printf("A      = %.5f +/- %.5f\n", position[0], chi*Math::sqrt(covar[0,0]/dof))
printf("lambda = %.5f +/- %.5f\n", position[1], chi*Math::sqrt(covar[1,1]/dof))
printf("b      = %.5f +/- %.5f\n", position[2], chi*Math::sqrt(covar[2,2]/dof))
