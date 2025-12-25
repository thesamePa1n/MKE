from dolfin import * 

mesh = Mesh('mesh.xml')
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml')
domains = MeshFunction('size_t', mesh, 'mesh_physical_region.xml')
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
dx = Measure('dx', domain=mesh, subdomain_data=domains)

f = Constant(1.0)
kinv = Constant(1.0)
g = Constant(0.0)
r_coeff = Constant(0.1)

V = FiniteElement('RT', mesh.ufl_cell(), 1)
P = FiniteElement('DG', mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, MixedElement(V, P))

(q, u) = TrialFunctions(W)
(r, v) = TestFunctions(W)

n = FacetNormal(mesh)

bc = DirichletBC(W.sub(1), g, boundaries, 2)

a = inner(kinv*q, r)*dx - div(r)*u*dx + div(q)*v*dx  + inner(q, n)*inner(r, n)*ds(2)
L = f * v * dx - g*inner(r, n)*ds(1) - Constant(0.1)*inner(r, n)*ds(2)

w = Function(W)
solve(a == L, w, bc) 

(q, u) = w.split()

fileq = File('/mnt/c/Users/Vasya/Downloads/q.pvd') 
fileu = File('/mnt/c/Users/Vasya/Downloads/u.pvd') 

fileq << q
fileu << u
