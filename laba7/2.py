from dolfin import * 

mesh = Mesh('mesh.xml')
boundaries = MeshFunction('size_t', mesh, 'mesh_facet_region.xml')

f = Constant(1.0)
kinv = Constant(1.0)
g = Constant(0.0)
r_coeff = Constant(0.1)

V = VectorFunctionSpace(mesh, 'CG', 2)
P = FunctionSpace(mesh, 'CG', 1)
W = MixedFunctionSpace(V, P)

(q, u) = TrialFunctions(W)
(r, v) = TestFunctions(W)

n = FacetNormal(mesh)

ds = Measure('ds')(domain=mesh, subdomain_data=boundaries)

bc = DirichletBC(W.sub_space(1), g, boundaries, 1)

a = inner(kinv*q, r)*dx + div(r)*u*dx + div(q)*v*dx + r_coeff*u*v*dx
L = f * v * dx + 0.1*inner(r, n)*ds(2)

w = Function(W)
solve(a == L, w, bc) 

(q, u) = w.split()

fileq = File('/mnt/c/Users/Vasya/Downloads/q.pvd') 
fileu = File('/mnt/c/Users/Vasya/Downloads/u.pvd') 

fileq << q
fileu << u
