from dolfin import * 

mesh = Mesh('mesh.xml') 
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml') 
boundaries = MeshFunction('size_t', mesh , 'mesh_facet_region.xml') 

V = FunctionSpace(mesh , 'Lagrange' , 1) 

dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

f = Constant(1.0)

bc_1 = DirichletBC(V, Constant(0.0), boundaries, 1)
bc_2 = DirichletBC(V, Constant(0.0), boundaries, 2)
bc_3 = DirichletBC(V, Constant(0.0), boundaries, 3)
bc_4 = DirichletBC(V, Constant(0.0), boundaries, 4)

bcs = [bc_1, bc_2, bc_3, bc_4]

k_1 = Constant(1) 
k_2 = Constant(2)
k_3 = Constant(820)
k_4 = Constant(420)

u = TrialFunction(V) 
v = TestFunction(V) 

a = k_1 * inner(grad(u), grad(v)) * dx(1) + k_2 * inner(grad(u), grad(v)) * dx(2) +\
  k_3 * inner(grad(u), grad(v)) * dx(3) + k_4 * inner(grad(u), grad(v)) * dx(4)
L = f * v * dx

u = Function(V)

solve(a == L, u, bcs) 

file = File('./poisson.pvd') 
file << u