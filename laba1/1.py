from dolfin import *

mesh = UnitSquareMesh(20, 20)
V = FunctionSpace(mesh, "CG", 1)

k = Constant(1.0)
r = Constant(0.1) 
f = Constant(1.0)
g_value = 0.0

u = TrialFunction(V)
v = TestFunction(V)

def g_boundary(x, on_boundary):
    return on_boundary and near(x[1], 0.0)

bc = DirichletBC(V, Constant(g_value), g_boundary)

a = k * inner(grad(u), grad(v)) * dx + r * u * v * dx + u * v * ds
L = f * v * dx + 0.1 * v * ds

u = Function(V)
solve(a == L, u, bc)

file = File("diffusion_reaction.pvd")
file << u
