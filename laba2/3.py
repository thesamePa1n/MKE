from dolfin import * 

mesh = Mesh('mesh.xml') 
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml') 
boundaries = MeshFunction ('size_t', mesh , 'mesh_facet_region.xml') 

V = FunctionSpace (mesh , 'Lagrange' , 1) 

dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

f = Constant(1.0)
g_L = Constant(70.0)
g_R = Constant(20.0)
g_B = Constant(-30.0)
g_T = Constant(10000.0)
alpha = Constant(150.0)
beta = Constant(150.0)

bc_1 = DirichletBC(V, g_L, boundaries, 1)
bc_2 = DirichletBC(V, g_R, boundaries, 2)

bcs = [bc_1, bc_2]

k_1 = Constant(1) 
k_2 = Constant(420)
k_3 = Constant(100)
k_4 = Constant(800)

u = TrialFunction(V) 
v = TestFunction(V) 

a = k_1 * inner(grad(u), grad(v)) * dx(1) + k_2 * inner(grad(u), grad(v)) * dx(2) + k_3 * inner(grad(u), grad(v)) * dx(3) + k_4 * inner(grad(u), grad(v)) * dx(4) + alpha * u * v * ds(3) + beta * u * v * ds(4)
L = f * v * dx(1) + f * v * dx(2) + alpha * g_B * v * ds(3) + beta * g_T * v * ds(4)

u = Function(V)

solve(a == L, u, bcs) 

file = File('./poisson2.pvd') 
file << u