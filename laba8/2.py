from dolfin import * 
import numpy

mesh = Mesh('mesh.xml') 
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml') 

def D(u):
	return 2 * u

def dD(u):
	return Constant(2.0)

f = Constant(0.0)
k = Expression("1 + pow(x[0], 2) + 2*pow(x[1], 2)", degree=2)
r = Constant(0.2)

V = FunctionSpace(mesh, 'Lagrange', 1) 
u = TrialFunction(V) 
v = TestFunction(V) 

bc1 = DirichletBC(V, Constant(10.0), boundaries, 1)
bc2 = DirichletBC(V, Constant(1.0), boundaries, 2)   
bcs = [bc1, bc2]

#Define variational problem for initial guess (q(u)=1, i.e., m=0)
a = k*inner(grad(u), grad(v))*dx + r*u*v*dx
L = f * v * dx
u_ = Function(V)
solve(a == L, u_, bcs)

bc1 = DirichletBC(V, Constant(0.0), boundaries, 1)
bc2 = DirichletBC(V, Constant(0.0), boundaries, 2)  
bcs = [bc1 , bc2]

# Define variational problem for Newton iteration
du = TrialFunction(V)
a = k*inner(grad(du), grad(v))*dx + r*du*v*dx + inner(D(u_)*grad(du), grad(v))*dx + inner(dD(u_)*du*grad(u_), grad(v))*dx
L = -k*inner(grad(u_), grad(v))*dx - r*u_*v*dx - inner(D(u_)*grad(u_), grad(v))*dx + f*v*dx

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