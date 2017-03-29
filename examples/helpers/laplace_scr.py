import numpy as np
from numba import autojit, jit
import scipy.io
import matplotlib.pyplot as plt
import copy
from scipy.sparse.linalg import gmres, LinearOperator,spilu
from scipy.sparse.linalg import spsolve
from time import time
from scipy.io import FortranFile
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree


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
    x, y, z = np.mgrid[0:N, 0:N, 0:N]
    data = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    tree = KDTree(data, leaf_size=2)
    prm = np.array(tree.idx_array)
    SE = csr_mat[:, prm][prm, :]
    return SE
def _laplace3d_4(d):
    N = (2 ** d)
    csr_mat = _csr_lap3d_4(N)
    x, y, z = np.mgrid[0:N, 0:N, 0:2*N]
    data = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    tree = KDTree(data, leaf_size=2)
    prm = np.array(tree.idx_array)
    SE = csr_mat[:, prm][prm, :]
    return SE
def _laplace3d_5(d):
    N = (2 ** d)
    csr_mat = _csr_lap3d_5(N)
    x, y, z = np.mgrid[0:N, 0:2*N, 0:2*N]
    data = np.array(zip(x.ravel(), y.ravel(), z.ravel()))
    tree = KDTree(data, leaf_size=2)
    prm = np.array(tree.idx_array)
    SE = csr_mat[:, prm][prm, :]
    return SE
def gen_3d_lap_good_prm(d):
    if d % 3 == 0:
        csr = _laplace3d(d/3)
    elif d % 3 == 1:
        csr = _laplace3d_4(d/3)
    elif d % 3 == 2:
        csr = _laplace3d_5(d/3)
    return csr
def _csr_lap2d(N):
    import  pyamg
    return scipy.sparse.csr_matrix(pyamg.gallery.laplacian.poisson((N,N)))
def _csr_lap2d_4(N):
    import  pyamg
    return scipy.sparse.csr_matrix(pyamg.gallery.laplacian.poisson((N,2*N)))
def _laplace2d(d):
    N = (2 ** d)
    csr_mat = _csr_lap2d(N)
    x, y = np.mgrid[0:N, 0:N]
    data = np.array(zip(x.ravel(), y.ravel()))
    tree = KDTree(data, leaf_size=2)
    prm = np.array(tree.idx_array)
    SE = csr_mat[:, prm][prm, :]
    return SE
def _laplace2d_4(d):
    N = (2 ** d)
    csr_mat = _csr_lap2d_4(N)
    x, y, = np.mgrid[0:N, 0:2*N]
    data = np.array(zip(x.ravel(), y.ravel()))
    tree = KDTree(data, leaf_size=2)
    prm = np.array(tree.idx_array)
    SE = csr_mat[:, prm][prm, :]
    return SE

def gen_2d_lap_good_prm(d):
    if d % 2 == 0:
        csr = _laplace3d(d/2)
    elif d % 2 == 1:
        csr = _laplace3d_4(d/2)
    return csr
