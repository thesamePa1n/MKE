from dolfin import *


mesh = UnitSquareMesh(3, 2) 

print("1. КООРДИНАТЫ СРЕДНИХ ТОЧЕК ЯЧЕЕК:")
for i in range(mesh.num_cells()):
    cell = Cell(mesh, i)
    mid = cell.midpoint()
    print(f"Ячейка {cell.index()}: ({mid.x():.3f}, {mid.y():.3f})")

print("\n2. КООРДИНАТЫ ВСЕХ УЗЛОВ СЕТКИ:")
for i in range(mesh.num_cells()): 
  cell = Cell(mesh, i) 
  vis = cell.entities(0) 
  for vi in vis: 
    vertex = Vertex(mesh, vi) 
    xx = vertex.point().x() 
    yy = vertex.point().y() 
    print(f"Узел({xx}, {yy})")

print("\n3. ИНДЕКСЫ УЗЛОВ ДЛЯ КАЖДОЙ ЯЧЕЙКИ:")
for cell in cells(mesh):
    vertices = cell.entities(0)  # 0 - вершины (узлы)
    print(f"Ячейка {cell.index()}: узлы {list(vertices)}")

# 4. Все узлы, лежащие на границе
print("\n4. УЗЛЫ НА ГРАНИЦЕ:")
boundary_nodes = set()
for facet in facets(mesh):
    if facet.exterior():
        vertices = facet.entities(0)
        boundary_nodes.update(vertices)

print(f"Граничные узлы: {sorted(boundary_nodes)}")

# Визуализация
print("\n=== ВИЗУАЛИЗАЦИЯ ===")
print("Сетка сохранена в файл structured_mesh.pvd")
File("structured_mesh.pvd") << mesh