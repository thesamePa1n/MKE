from dolfin import * 

N = 32
mesh = UnitSquareMesh(N, N)

k = Constant(1.0)
f = Constant(1.0)
g = Constant(0.0)
r = Constant(0.1)

class DirichletBoundary(SubDomain): 
  def inside(self, x, on_boundary):
	  return on_boundary and (x[0] > (1.0 - DOLFIN_EPS) or x[1] > (1.0 - DOLFIN_EPS))

class RobinBoundary(SubDomain): 
	def inside(self, x, on_boundary): 
		return on_boundary and (x[0] < DOLFIN_EPS or x[1] < DOLFIN_EPS)

V = FunctionSpace(mesh, 'Lagrange', 1) 
boundaries = MeshFunction("size_t", mesh, 1)
RobinBoundary().mark(boundaries, 2) 
DirichletBoundary().mark(boundaries, 1) 

ds = Measure('ds')(domain=mesh, subdomain_data=boundaries)

u = TrialFunction(V) 
v = TestFunction(V) 

n = FacetNormal(mesh)
h = CellDiameter(mesh)
h_avg = (h('+')+h('-'))/2
alpha = 4.0
gamma = 8.0

a = inner(k*grad(u), grad(v)) * dx + r*u*v*dx\
  - inner(avg(k*grad(v)), jump(u, n)) * dS \
  - inner(jump(v, n), avg(k*grad(u))) * dS \
  + avg(k)*alpha/h_avg*inner(jump(v, n), jump(u, n)) * dS \
  - inner(k*grad(v), u*n)*ds(1) \
  - inner(v*n, k*grad(u))*ds(1) \
  + k*(gamma/h)*u*v*ds(1) + u*v*ds(2)
L = f*v*dx - g*inner(k*grad(v), n)*ds(1) + k*(gamma/h)*g*v*ds(1) + 0.1*v*ds(2)

u = Function(V)
solve(a == L, u) 

file = File('/mnt/c/Users/Vasya/Downloads/laba7-1.pvd') 
file << u
