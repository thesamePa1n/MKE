from dolfin import *
import time

start_time = time.time()

mesh = Mesh('mesh.xml') 
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml') 
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml') 

V = FunctionSpace(mesh, 'Lagrange', 1) 

dx = Measure('dx', domain=mesh, subdomain_data=domains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

f = Constant(1.0)
k = Constant(7.0)

bc_1 = DirichletBC(V, Constant(0.0), boundaries, 1)
bc_2 = DirichletBC(V, Constant(0.0), boundaries, 2)
bcs = [bc_1, bc_2]

u = TrialFunction(V) 
v = TestFunction(V)

a = k * inner(grad(u), grad(v)) * dx
L = f * v * dx

u = Function(V)
U = u.vector()
A = assemble(a)
b = assemble(L)

for bc in bcs:
  bc.apply(A, b)
  
solver = KrylovSolver('bicgstab', 'hypre_amg')

cg_prm = solver.parameters
cg_prm["absolute_tolerance"] = 1E-7     
cg_prm["relative_tolerance"] = 1E-4       
cg_prm["maximum_iterations"] = 1000 

iter = solver.solve(A, U, b)

end_time = time.time()

print(f'iter = {iter}')
print(f"Время выполнения: {end_time-start_time:.3f} секунд")

File('/mnt/c/Users/Vasya/Downloads/laba6.pvd') << u