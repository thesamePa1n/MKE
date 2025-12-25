from dolfin import * 
import numpy

mesh = Mesh('mesh.xml') 
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml') 
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml')
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx = Measure('dx', domain=mesh, subdomain_data=domains)

def k(u):
	return exp(u)

f = Constant(0.0)

V = FunctionSpace(mesh, 'Lagrange', 1) 
u = TrialFunction(V) 
v = TestFunction(V) 

bc1 = DirichletBC(V, Constant(5.0), boundaries, 1)
bcs = [bc1]

uk = interpolate(Expression(("0.0"), degree=1), V)
a = inner(k(uk)*grad(u), grad(v))*dx
L = f * v * dx + 2*v*ds(2)

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
