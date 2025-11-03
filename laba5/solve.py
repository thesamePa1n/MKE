from dolfin import * 

mesh1 = Mesh('mesh1.xml') 
domains1 = MeshFunction('size_t', mesh1, 'mesh1_physical_region.xml') 
boundaries1 = MeshFunction('size_t', mesh1, 'mesh1_facet_region.xml') 

mesh2 = Mesh('mesh2.xml') 
domains2 = MeshFunction('size_t', mesh2, 'mesh2_physical_region.xml') 
boundaries2 = MeshFunction('size_t', mesh2, 'mesh2_facet_region.xml') 

mesh3 = Mesh('mesh3.xml') 
domains3 = MeshFunction('size_t', mesh3, 'mesh3_physical_region.xml') 
boundaries3 = MeshFunction('size_t', mesh3, 'mesh3_facet_region.xml') 

mesh4 = Mesh('mesh4.xml') 
domains4 = MeshFunction('size_t', mesh4, 'mesh4_physical_region.xml') 
boundaries4 = MeshFunction('size_t', mesh4, 'mesh4_facet_region.xml') 

V1 = FunctionSpace(mesh1, 'Lagrange', 1) 
V2 = FunctionSpace(mesh2, 'Lagrange', 1) 
V3 = FunctionSpace(mesh3, 'Lagrange', 1) 
V4 = FunctionSpace(mesh4, 'Lagrange', 1) 

dx1 = Measure('dx', domain=mesh1, subdomain_data=domains1)
dx2 = Measure('dx', domain=mesh2, subdomain_data=domains2)
dx3 = Measure('dx', domain=mesh3, subdomain_data=domains3)
dx4 = Measure('dx', domain=mesh4, subdomain_data=domains4)

f = Constant(0.0)
g_L = Constant(10.0)
g_B = Constant(0.0)
g_T = Constant(0.0)
k = Constant(0.1)

bc_1_1 = DirichletBC(V1, g_L, boundaries1, 1)
bc_2_1 = DirichletBC(V1, g_B, boundaries1, 3)
bc_3_1 = DirichletBC(V1, g_T, boundaries1, 4)
bcs1 = [bc_1_1, bc_2_1, bc_3_1]

bc_1_2 = DirichletBC(V2, g_L, boundaries2, 1)
bc_2_2 = DirichletBC(V2, g_B, boundaries2, 3)
bc_3_2 = DirichletBC(V2, g_T, boundaries2, 4)
bcs2 = [bc_1_2, bc_2_2, bc_3_2]

bc_1_3 = DirichletBC(V3, g_L, boundaries3, 1)
bc_2_3 = DirichletBC(V3, g_B, boundaries3, 3)
bc_3_3 = DirichletBC(V3, g_T, boundaries3, 4)
bcs3 = [bc_1_3, bc_2_3, bc_3_3]

bc_1_4 = DirichletBC(V4, g_L, boundaries4, 1)
bc_2_4 = DirichletBC(V4, g_B, boundaries4, 3)
bc_3_4 = DirichletBC(V4, g_T, boundaries4, 4)
bcs4 = [bc_1_4, bc_2_4, bc_3_4]

u0_val = Constant(5.0)

u0_1 = interpolate(u0_val, V1)
u0_2 = interpolate(u0_val, V2)
u0_3 = interpolate(u0_val, V3)
u0_4 = interpolate(u0_val, V4)

T = 6.0
tau = 0.06

u1 = TrialFunction(V1) 
v1 = TestFunction(V1) 

u2 = TrialFunction(V2) 
v2 = TestFunction(V2) 

u3 = TrialFunction(V3) 
v3 = TestFunction(V3) 

u4 = TrialFunction(V4) 
v4 = TestFunction(V4) 

a1 = (u1/tau) * v1 * dx1 + k * inner(grad(u1), grad(v1)) * dx1
L1 = (u0_1/tau) * v1 * dx1

a2 = (u2/tau) * v2 * dx2 + k * inner(grad(u2), grad(v2)) * dx2
L2 = (u0_2/tau) * v2 * dx2

a3 = (u3/tau) * v3 * dx3 + k * inner(grad(u3), grad(v3)) * dx3
L3 = (u0_3/tau) * v3 * dx3

a4 = (u4/tau) * v4 * dx4 + k * inner(grad(u4), grad(v4)) * dx4
L4 = (u0_4/tau) * v4 * dx4

u1 = Function(V1)
u2 = Function(V2)
u3 = Function(V3)
u4 = Function(V4)

file1 = File('/mnt/c/Users/Vasya/Downloads/time_dep1.pvd')
file2 = File('/mnt/c/Users/Vasya/Downloads/time_dep2.pvd')
file3 = File('/mnt/c/Users/Vasya/Downloads/time_dep3.pvd')
file4 = File('/mnt/c/Users/Vasya/Downloads/time_dep4.pvd')

t = 0
cnt = 0
while t < T:
    t += tau
    solve(a1 == L1, u1, bcs1) 
    solve(a2 == L2, u2, bcs2) 
    solve(a3 == L3, u3, bcs3) 
    solve(a4 == L4, u4, bcs4) 
    file1 << u1 
    file2 << u2
    file3 << u3
    file4 << u4
    u0_1.assign(u1)
    u0_2.assign(u2)
    u0_3.assign(u3)
    u0_4.assign(u4)

    if (cnt % 10 == 0 and cnt != 0):
        u_e1 = interpolate(u4, V1)
        u_e2 = interpolate(u4, V2)
        u_e3 = interpolate(u4, V3)
        
        E1 = (u_e1 - u1)*dx1
        E2 = (u_e2 - u2)*dx2
        E3 = (u_e3 - u3)*dx3
        E1_error1 = abs(assemble(E1))*100
        E1_error2 = abs(assemble(E2))*100
        E1_error3 = abs(assemble(E3))*100
        print(f"abs_error: E1_error1={E1_error1}%, E1_error2={E1_error2}%, E1_error3={E1_error3}%")
        
        EL2_t1 = inner(u1 - u_e1, u1 - u_e1) * dx1
        EL2_t2 = inner(u2 - u_e2, u2 - u_e2) * dx2
        EL2_t3 = inner(u3 - u_e3, u3 - u_e3) * dx3
        EL2_b1 = inner(u1, u1) * dx1
        EL2_b2 = inner(u2, u2) * dx2
        EL2_b3 = inner(u3, u3) * dx3
        EL2_error1 = sqrt(abs(assemble(EL2_t1)) / abs(assemble(EL2_b1)))*100    
        EL2_error2 = sqrt(abs(assemble(EL2_t2)) / abs(assemble(EL2_b2)))*100    
        EL2_error3 = sqrt(abs(assemble(EL2_t3)) / abs(assemble(EL2_b3)))*100
        print(f"error L2: EL2_error1={EL2_error1}%, EL2_error2={EL2_error2}%, EL2_error3={EL2_error3}%")

        EH1_t1 = inner(grad(u1 - u_e1), grad(u1 - u_e1)) * dx1
        EH1_t2 = inner(grad(u2 - u_e2), grad(u2 - u_e2)) * dx2
        EH1_t3 = inner(grad(u3 - u_e3), grad(u3 - u_e3)) * dx3
        EH1_b1 = inner(grad(u1), grad(u1)) * dx1
        EH1_b2 = inner(grad(u2), grad(u2)) * dx2
        EH1_b3 = inner(grad(u3), grad(u3)) * dx3
        EH1_error1 = sqrt(abs(assemble(EH1_t1)) / abs(assemble(EH1_b1)))*100    
        EH1_error2 = sqrt(abs(assemble(EH1_t2)) / abs(assemble(EH1_b2)))*100    
        EH1_error3 = sqrt(abs(assemble(EH1_t3)) / abs(assemble(EH1_b3)))*100
        print(f"error H1: EH1_error1={EH1_error1}%, EH1_error2={EH1_error2}%, EH1_error3={EH1_error3}%")   
    cnt += 1


