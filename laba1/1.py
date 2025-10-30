from dolfin import *

mesh = UnitSquareMesh(6, 4)
V = FunctionSpace(mesh, "CG", 1)
g = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)

def g_boundary(x, on_boundary):
  return on_boundary

bc = DirichletBC(V, g, g_boundary)
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(grad(u), grad(v)) * dx
L = f * v * dx
u = Function(V)
solve(a == L, u, bc)
file = File("mesh.pvd")
file << mesh
file2 = File("poisson.pvd")
file2 << u