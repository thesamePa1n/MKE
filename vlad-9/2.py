from dolfin import *
import numpy

parameters["reorder_dofs_serial"] = False

N = 5
mesh = UnitIntervalMesh(N)
k = 0.1
f = 0
Nc = mesh.num_cells()
NN = mesh.num_vertices()
AA = numpy.zeros((NN, NN)) 
FF = numpy.zeros(NN)

for ci in range(Nc):
  cell = Cell(mesh, ci)
  vertices = cell.entities(0)
  vert0 = Vertex(mesh, vertices[0])
  vert1 = Vertex(mesh, vertices[1])
  h = abs(vert0.point().x() - vert1.point().x())
  dof0 = vertices[0]
  dof1 = vertices[1]
  AA[dof0][dof0] += 1.0 / h * k
  AA[dof1][dof0] += -1.0/ h * k
  AA[dof0][dof1] += -1.0/ h * k
  AA[dof1][dof1] += 1.0 / h * k
  
  FF[dof0] += f * h / 2
  FF[dof1] += f * h / 2
  
AA[0, :] = 0
AA[-1, :] = 0
AA[0, 0] = 1.0
AA[-1, -1] = 1.0
FF[0] = 5.0
FF[-1] = 0.0
  
U = numpy.linalg.solve(AA, FF)
print("Our solve")
print(U)
print()

V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
kk = Constant(k)
ff = Constant(f)
a = inner(kk * grad(u), grad(v)) * dx
L = ff * v * dx

def left_boundary(x, on_boundary):
    return on_boundary and x[0] < DOLFIN_EPS

def right_boundary(x, on_boundary):
    return on_boundary and x[0] > 1.0 - DOLFIN_EPS

bc_left = DirichletBC(V, Constant(5.0), left_boundary)
bc_right = DirichletBC(V, Constant(0.0), right_boundary)

A = Matrix()
b = Vector()
assemble(a, A)
assemble(L, b)
bc_left.apply(A, b)
bc_right.apply(A, b)

u = Function(V)
solve(A, u.vector(), b)

print('Fenics solve')
print(u.vector().get_local())