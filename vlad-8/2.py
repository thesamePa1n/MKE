from dolfin import * 
import numpy

mesh = Mesh('mesh.xml') 
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml') 
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml')
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx = Measure('dx', domain=mesh, subdomain_data=domains)

def k(u):
	return exp(u)

def dK(u):
	return exp(u)

f = Constant(0.0)

V = FunctionSpace(mesh, 'Lagrange', 1) 
u = TrialFunction(V) 
v = TestFunction(V) 

bc1 = DirichletBC(V, Constant(5.0), boundaries, 1)
bcs = [bc1]

#Define variational problem for initial guess (k(u)=1, i.e., m=0)
a = inner(grad(u), grad(v))*dx
L = f * v * dx + 2*v*ds(2)
u_ = Function(V)
solve(a == L, u_, bcs)

bc1 = DirichletBC(V, Constant(0.0), boundaries, 1)
bcs = [bc1]

# Define variational problem for Newton iteration
du = TrialFunction(V)
a = inner(k(u_)*grad(du), grad(v))*dx + inner(dK(u_)*du*grad(u_), grad(v))*dx
L = -inner(k(u_)*grad(u_), grad(v))*dx + f*v*dx + 2*v*ds(2)

#Newton iteration at the PDE level
du = Function(V)
u = Function(V)
omega = 1
eps = 1.0
tol = 1.0E-5
iter = 0
maxiter = 25

while eps > tol and iter < maxiter:
	iter += 1
	print(iter, 'iteration', end=' ')
	solve(a == L, du, bcs) 
	eps = numpy.linalg.norm(du.vector().get_local(), ord=numpy.Inf)
	print('Norm:', eps)
	u.vector()[:] = u_.vector() + omega * du.vector()
	u_.assign(u)

file = File('/mnt/c/Users/Vasya/Downloads/laba8-2.pvd') 
file << u