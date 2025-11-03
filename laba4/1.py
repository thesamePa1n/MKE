from dolfin import *


mesh = UnitSquareMesh(3, 2) 
V = FunctionSpace(mesh, "CG", 1)

print("\nкоординаты всех средних точек ячеек")
for i in range(mesh.num_cells()):
    cell = Cell(mesh, i)
    mid = cell.midpoint()
    print(f"Ячейка {cell.index()}: ({mid.x():.3f}, {mid.y():.3f})")

print("\nкоординаты узлов сетки")
for i in range(mesh.num_cells()): 
    cell = Cell(mesh, i) 
    vis = cell.entities(0) 
    for vi in vis: 
        vertex = Vertex(mesh, vi) 
        xx = vertex.point().x() 
        yy = vertex.point().y() 
        print(f"узел({xx}, {yy})")

print("\nдля каждой ячейки индексы узлов")
for i in range(mesh.num_cells()):
    cell = Cell(mesh, i)
    vis = cell.entities(0)
    print(f"Ячейка {cell.index()}: узлы {list(vis)}")

print("\nвсех узлов, лежащих на границе")
for vi in range(mesh.num_vertices()):
    vertex = Vertex(mesh, vi)
    xx = vertex.point().x()
    yy = vertex.point().y()
    
    if xx < DOLFIN_EPS:
        print(f'левая граница: узел({xx}, {yy})')
    elif abs(yy - 1.0) < DOLFIN_EPS:
        print(f'верхняя граница: узел({xx}, {yy})')
    elif abs(xx - 1.0) < DOLFIN_EPS:
        print(f'правая граница: узел({xx}, {yy})')
    elif yy < DOLFIN_EPS:
        print(f'нижняя граница: узел({xx}, {yy})')

File("structured_mesh.pvd") << mesh

u = TrialFunction(V)
v = TestFunction(V)

k = Constant(0.3)
f = Constant(1.0)
g = Constant(0.0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, g, boundary)

a = k * inner(grad(u), grad(v)) * dx
L = f * v * dx

u = Function(V)

solve(a == L, u, bc)
File("solution.pvd") << u

vert_to_dof = vertex_to_dof_map(V)
u_array = u.vector().get_local()
print(f"общее количество DOF: {len(u_array)}")
print()

for vi in range(mesh.num_vertices()):
    dof = vert_to_dof[vi]
    vertex = Vertex(mesh, vi)
    xx = vertex.point().x()
    yy = vertex.point().y()
    print(f"{dof}, ({xx:.2f}, {yy}), {u_array[dof]}")
print()

cell_idx = 2
cell = Cell(mesh, cell_idx)
vertices = cell.entities(0)
dofs = [vert_to_dof[vi] for vi in vertices]
values = [u_array[dof] for dof in dofs]

print(f"Ячейка {cell_idx}: DOF {dofs}, значения {values}")

average = sum(values) / len(values)
print("Среднее значение", average)