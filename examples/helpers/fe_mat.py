import numpy as np
import scipy as sp
from fenics import *
from mshr import *
#import  mshr


def complex_mesh_2d(r=100):
    # Create circles as Circle(Center, Radius)
    circle1 = Circle(Point(0,0), 5)
    circle2 = Circle(Point(-1,0), 1)


    domain = circle1# - circle2
    mesh = generate_mesh(domain, r)

    #plot(mesh, interactive=True)
    #plt.show()
    V = FunctionSpace(mesh, 'P', 1)
    # Define boundary condition
    u_D = Expression("x[0] - x[1]", degree=2)
    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, u_D, boundary)
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v))*dx
    L = f*v*dx
    A, b = assemble_system(a, L, bc)
    A_mat = as_backend_type(A).mat()
    from scipy.sparse import csr_matrix
    A_sparray = csr_matrix(A_mat.getValuesCSR()[::-1], shape = A_mat.size)
    #print A_sparray.shape
    #plt.spy(A_sparray,markersize=0.1)
    rhs = b.array()
    return mesh, A_sparray, rhs



def complex_mesh_3d(r = 60):
    sp1 = Sphere(Point(0,0,0),1)
    sp2 = Sphere(Point(0,0,0),0.5)
    dom = sp1 - sp2
    mesh = generate_mesh(dom, r)
    parameters['reorder_dofs_serial'] = False
    #plot(mesh, interactive=True)
    #plt.show()
    V = FunctionSpace(mesh, 'P', 1)
    # Define boundary condition
    u_D = Expression("x[0] - x[1]", degree=2)
    def boundary(x, on_boundary):
        return on_boundary
    bc = DirichletBC(V, u_D, boundary)
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    f = Constant(-6.0)
    a = dot(grad(u), grad(v))*dx
    L = f*v*dx
    A, b = assemble_system(a, L, bc)
    A_mat = as_backend_type(A).mat()
    from scipy.sparse import csr_matrix
    A_sparray = csr_matrix(A_mat.getValuesCSR()[::-1], shape = A_mat.size)
    #print A_sparray.shape
    #plt.spy(A_sparray,markersize=0.1)
    rhs = b.array()
    return mesh, A_sparray, rhs