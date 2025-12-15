from dolfin import *
import numpy

parameters["reorder_dofs_serial"] = False

N = 5
mesh = UnitIntervalMesh(N)
k = 0.7
q = 9.0
Nc = mesh.num_cells()
NN = mesh.num_vertices()
AA = numpy.zeros((NN, NN)) 

for ci in range(Nc):
  cell = Cell(mesh, ci)
  vertices = cell.entities(0)
  vert0 = Vertex(mesh, vertices[0])
  vert1 = Vertex(mesh, vertices[1])
  h = abs(vert0.point().x() - vert1.point().x())
  dof0 = vertices[0]
  dof1 = vertices[1]
  AA[dof0][dof0] += 1.0 / h * k + h / 6 * 2 * q
  AA[dof1][dof0] += -1.0/ h * k + h/6*q
  AA[dof0][dof1] += -1.0/ h * k + h/6*q
  AA[dof1][dof1] += 1.0 / h * k + h / 6 * 2 * q
  
print(AA)
print() 

V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)
kk = Constant(k)
qq = Constant(q)
a = inner(kk * grad(u), grad(v)) * dx + qq * u * v * dx
A = Matrix()
assemble(a, A)
print(A.array()) 
