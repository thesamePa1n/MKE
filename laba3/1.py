from dolfin import * 

mesh = Mesh('mesh.xml') 
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml') 
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml') 

V = FunctionSpace(mesh, 'Lagrange', 1) 

dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

f = Constant(0.0)
g_L = Constant(-10.0)
g_R = Constant(15.0)
g_B = Constant(0.0)
g_T = Constant(1.0)
alpha = Constant(100.0)

bc_1 = DirichletBC(V, g_L, boundaries, 1)
bc_2 = DirichletBC(V, g_R, boundaries, 2)
bcs = [bc_1, bc_2]

u0_val = Constant(5.0)
u0 = interpolate(u0_val, V)

k_1 = Constant(50) 
k_2 = Constant(50)
k_3 = Constant(200)
k_4 = Constant(200)

C_1 = Constant(150)
C_2 = Constant(150)
C_3 = Constant(200)
C_4 = Constant(200)

T = 600.0
N = 10
tau = T / N

u = TrialFunction(V) 
v = TestFunction(V) 

a = (C_1/tau)*u*v*dx(1) + (C_2/tau)*u*v*dx(2) + (C_3/tau)*u*v*dx(3) + (C_4/tau)*u*v*dx(4) +\
    k_1 * inner(grad(u), grad(v)) * dx(1) + k_2 * inner(grad(u), grad(v)) * dx(2) +\
    k_3 * inner(grad(u), grad(v)) * dx(3) + k_4 * inner(grad(u), grad(v)) * dx(4) +\
    alpha*u*v*ds(4) + alpha*u*v*ds(3)
L = (C_1/tau)*u0*v*dx(1) + (C_2/tau)*u0*v*dx(2) + (C_3/tau)*u0*v*dx(3) + (C_4/tau)*u0*v*dx(4) +\
    alpha*g_B*v*ds(4) + alpha*g_T*v*ds(3)

u = Function(V)
file = File('/mnt/c/Users/Vasya/Downloads/time_dep.pvd')

t = 0
while t < T:
	t += tau / 10
	solve(a == L, u, bcs) 
	file << u 
	u0.assign(u)
