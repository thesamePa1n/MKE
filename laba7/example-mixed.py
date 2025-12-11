from dolfin import * 

mesh = Mesh('mesh2.xml')
boundaries = MeshFunction('size_t', mesh, 'mesh2_facet_region.xml')

f = Constant(1.0)
kinv = Constant(1.0)
g2 = Constant((0.0, 0.0))
g1 = Constant(0.0)

# V = VectorFunctionSpace(mesh, 'CG', 2)
# P = FunctionSpace(mesh, 'CG', 1)
# W = V*P

V = VectorElement("CG", mesh.ufl_cell(), 2) 
P = FiniteElement("CG", mesh.ufl_cell(), 1) 
WE = MixedElement([V, P])                    
W = FunctionSpace(mesh, WE)  

(q, u) = TrialFunctions(W)
(r, v) = TestFunctions(W)

n = FacetNormal(mesh)

ds = Measure('ds')(domain=mesh, subdomain_data=boundaries)

bc = DirichletBC(W.sub(0), g2, boundaries, 2)

a = inner(kinv*q, r)*dx + div(r)*u*dx + div(q)*v*dx
L = f * v * dx + g1*inner(r, n)*ds(1)

w = Function(W)
solve(a == L, w, bc) 

(q, u) = w.split()

fileq = File('/mnt/c/Users/Vasya/Downloads/q.pvd') 
fileu = File('/mnt/c/Users/Vasya/Downloads/u.pvd') 

fileq << q
fileu << u
