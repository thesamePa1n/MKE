from dolfin import * 
import numpy, sys


N = 32
mesh = UnitSquareMesh(N, N)
m = 2
def q(u):
	return (1 + u)**m

f = Constant(0.0)
g1 = Constant(1.0)
g2 = Constant(0.0)

def g1_boundary(x): 
	return x[0] < DOLFIN_EPS

def g2_boundary(x): 
	return x[0] > 1.0 - DOLFIN_EPS

V = FunctionSpace(mesh, 'Lagrange', 1) 
u = TrialFunction(V) 
v = TestFunction(V) 

bc1 = DirichletBC(V, g1, g1_boundary) 
bc2 = DirichletBC(V, g2, g2_boundary)   
bcs = [bc1, bc2]
uk = interpolate(Expression(("0.0"), degree=1), V)
a = inner(q(uk)*grad(u), grad(v)) * dx 
L = f * v * dx

#Picard iterations
u = Function(V)
eps = 1.0
tol = 1.0E-5
iter = 0
maxiter = 25

while eps > tol and iter < maxiter:
	iter += 1
	solve(a == L, u, bcs)
	# u = problem.solve()
	diff = u.vector().get_local() - uk.vector().get_local()
	eps = numpy.linalg.norm(diff, ord=numpy.Inf)
	print('Norm, iter=%d: %g' % (iter, eps))
	uk.assign(u)

file = File('/mnt/c/Users/Vasya/Downloads/poisson.pvd') 
file << u
# Find max error
u_exact = Expression("pow((pow(2, m+1)-1)*x[0] + 1, 1.0/(m+1)) - 1", m=m, degree=1)
u_e = interpolate(u_exact, V)
diff = numpy.abs(u_e.vector().get_local() - u.vector().get_local()).max()
print('Max error:', diff)
