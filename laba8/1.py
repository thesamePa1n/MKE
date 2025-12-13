from dolfin import * 
import numpy

mesh = Mesh('mesh.xml') 
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml') 

def D(u):
	return 2 * u

f = Constant(0.0)
k = Expression("1 + pow(x[0], 2) + 2*pow(x[1], 2)", degree=2)
r = Constant(0.2)

V = FunctionSpace(mesh, 'Lagrange', 1) 
u = TrialFunction(V) 
v = TestFunction(V) 

bc1 = DirichletBC(V, Constant(10.0), boundaries, 1)
bc2 = DirichletBC(V, Constant(1.0), boundaries, 2)  
bcs = [bc1, bc2]

uk = interpolate(Expression(("0.0"), degree=1), V)
a = k*inner(grad(u), grad(v))*dx + D(uk)*inner(grad(u), grad(v))*dx + r*u*v*dx
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

file = File('/mnt/c/Users/Vasya/Downloads/laba8-1.pvd') 
file << u
