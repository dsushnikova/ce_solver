import numpy as np
from numba import autojit, jit
import scipy.io
#%matplotlib inline
import matplotlib.pyplot as plt
import copy
from scipy.sparse.linalg import gmres, LinearOperator,spilu
from scipy.sparse.linalg import spsolve
from time import time
from scipy.io import FortranFile
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
#from scikits.sparse.cholmod import cholesky

#import sys
#sys.path.append("~/INM/fdirsparsefmm/py_ce/fortran_core/")

def _csr_lap3d(N):
    import  pyamg
    return scipy.sparse.csr_matrix(pyamg.gallery.laplacian.poisson( (N,N,N) ))
def _csr_lap3d_4(N):
    import  pyamg
    return scipy.sparse.csr_matrix(pyamg.gallery.laplacian.poisson( (N,N,2*N) ))
def _csr_lap3d_5(N):
    import  pyamg
    return scipy.sparse.csr_matrix(pyamg.gallery.laplacian.poisson( (N,2*N,2*N) ))
def _laplace3d(d):
    N = (2 ** d)
    csr_mat = _csr_lap3d(N)
    #x, y, z = np.mgrid[0:N, 0:N, 0:N]
    #data = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    #tree = KDTree(data, leaf_size=2)
    #prm = np.array(tree.idx_array)
    #SE = csr_mat[:, prm][prm, :]
    return _csr_lap3d(N)
def _laplace3d_4(d):
    N = (2 ** d)
    #sr_mat = _csr_lap3d_4(N)
    #x, y, z = np.mgrid[0:N, 0:N, 0:2*N]
    #ata = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    #tree = KDTree(data, leaf_size=2)
    #prm = np.array(tree.idx_array)
    #SE = csr_mat[:, prm][prm, :]
    return _csr_lap3d_4(N)
def _laplace3d_5(d):
    N = (2 ** d)
    #csr_mat = _csr_lap3d_5(N)
    #x, y, z = np.mgrid[0:N, 0:2*N, 0:2*N]
    #data = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    #tree = KDTree(data, leaf_size=2)
    #prm = np.array(tree.idx_array)
    #SE = csr_mat[:, prm][prm, :]
    return _csr_lap3d_5(N)

def gen_3d_lap(d):
    if d % 3 == 0:
        csr = _laplace3d(d/3)
    elif d % 3 == 1:
        csr = _laplace3d_4(d/3)
    elif d % 3 == 2:
        csr = _laplace3d_5(d/3)
    return csr

def make_coord(d):
    N = (2 ** (d/3))
    #print N
    if d % 3 == 0:
        x, y, z = np.mgrid[0:N, 0:N, 0:N]
    if d % 3 == 1:
        x, y, z = np.mgrid[0:N, 0:N, 0:2*N]
        #print x.shape
    if d % 3 == 2:
        x, y, z = np.mgrid[0:N, 0:2*N, 0:2*N]
    data = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    #print data.shape
    return data

def make_edges(d):
    N = (2 ** (d/3))
    if d % 3 == 0:
        ed = []
        for i in xrange(N**3):
            ed.append([])
        
        for x in xrange(1,N):
            for y in xrange(0,N):
                for z in xrange(0,N):
                    ed[x*N*N + N*y +z].append((x-1)*N*N+ N*y + z)
        for x in xrange(0,N):
            for y in xrange(1,N):
                for z in xrange(0,N):
                    ed[x*N*N + N*y +z].append(x*N*N+ N*(y-1) + z)
        for x in xrange(0,N):
            for y in xrange(0,N):
                for z in xrange(1,N):
                    ed[x*N*N + N*y +z].append(x*N*N+ N*y + z-1)
    
        for x in xrange(0,N-1):
            for y in xrange(0,N):
                for z in xrange(0,N):
                    ed[x*N*N + N*y +z].append((x+1)*N*N+ N*y + z)
        for x in xrange(0,N):
            for y in xrange(0,N-1):
                for z in xrange(0,N):
                    ed[x*N*N + N*y +z].append(x*N*N+ N*(y+1) + z)
        for x in xrange(0,N):
            for y in xrange(0,N):
                for z in xrange(0,N-1):
                    ed[x*N*N + N*y +z].append(x*N*N+ N*y + z+1)
    if d % 3 == 1:
        ed = []
        for i in xrange((N**3)*2):
            ed.append([])
        
        for x in xrange(1,N):
            for y in xrange(0,N):
                for z in xrange(0,N*2):
                    #print x,y,z
                    #print x*N*N*2 + N*2*y +z, (x-1)*N*N*2+ N*2*y + z
                    ed[x*N*N*2 + N*2*y +z].append((x-1)*N*N*2+ N*2*y + z)
        for x in xrange(0,N):
            for y in xrange(1,N):
                for z in xrange(0,N*2):
                    ed[x*N*N*2 + N*y*2 +z].append(x*N*N*2+ N*(y-1)*2 + z)
        for x in xrange(0,N):
            for y in xrange(0,N):
                for z in xrange(1,N*2):
                    ed[x*N*N*2 + N*2*y +z].append(x*2*N*N+ N*2*y + z-1)
    
        for x in xrange(0,N-1):
            for y in xrange(0,N):
                for z in xrange(0,N*2):
                    ed[x*N*N*2 + N*y*2 +z].append((x+1)*2*N*N+ N*2*y + z)
        for x in xrange(0,N):
            for y in xrange(0,N-1):
                for z in xrange(0,N*2):
                    ed[x*N*N*2 + N*y*2 +z].append(x*N*N*2+ N*(y+1)*2 + z)
        for x in xrange(0,N):
            for y in xrange(0,N):
                for z in xrange(0,N*2-1):
                    ed[x*N*N*2 + N*y*2 +z].append(x*N*N*2 + N*y*2 + z+1)
    if d % 3 == 2:
        ed = []
        for i in xrange((N**3)*4):
            ed.append([])
        
        for x in xrange(1,N):
            for y in xrange(0,N*2):
                for z in xrange(0,N*2):
                    #print x,y,z
                    #print x*N*N*4 + N*2*y +z, (x-1)*N*N*4+ N*2*y + z
                    ed[x*N*N*4 + N*2*y +z].append((x-1)*N*N*4+ N*2*y + z)
        for x in xrange(0,N):
            for y in xrange(1,N*2):
                for z in xrange(0,N*2):
                    ed[x*N*N*4 + N*y*2 +z].append(x*N*N*4+ N*(y-1)*2 + z)
        for x in xrange(0,N):
            for y in xrange(0,N*2):
                for z in xrange(1,N*2):
                    ed[x*N*N*4 + N*2*y +z].append(x*4*N*N+ N*2*y + z-1)
    
        for x in xrange(0,N-1):
            for y in xrange(0,N*2):
                for z in xrange(0,N*2):
                    ed[x*N*N*4 + N*y*2 +z].append((x+1)*4*N*N+ N*2*y + z)
        for x in xrange(0,N):
            for y in xrange(0,N*2-1):
                for z in xrange(0,N*2):
                    ed[x*N*N*4 + N*y*2 +z].append(x*N*N*4+ N*(y+1)*2 + z)
        for x in xrange(0,N):
            for y in xrange(0,N*2):
                for z in xrange(0,N*2-1):
                    ed[x*N*N*4 + N*y*2 +z].append(x*N*N*4 + N*y*2 + z+1)

    return ed

