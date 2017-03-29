import numpy as np
from numba import autojit, jit
import scipy.io
import matplotlib.pyplot as plt
import copy
from scipy.sparse.linalg import gmres, LinearOperator, spilu, cg
from scipy.sparse.linalg import spsolve
from time import time
from scipy.io import FortranFile
import matplotlib.pyplot as plt
import fortran_core
from fortran_core.core import prec
import inspect


def h2_dir(csr, rhs, eps=1e-3, r=3, block_size=8, bl_ar=[0],level=6, join=2, mj=[[0]]):
    if bl_ar[0] == 0:
        M = csr.shape[0]/block_size
        bl_ar = np.ones(M, order='F')*8
    else:
        M = bl_ar.shape[0]
    t0 = time()

    if mj[0][0] == 0:
        M_j = []
        M_j.append(M)
        mj = []
        for i in xrange(level-1):
            if M_j[-1]%join == 0:
                M_j.append(M_j[-1]/join)
                for i in xrange(M_j[-1]):
                    mj.append(join)
            else:
                M_j.append(M_j[-1]/join+1)
                for i in xrange(M_j[-1] - 1):
                    mj.append(join)
                mj.append(M_j[-2]%join)
        hyper_join = np.array(mj, order='F')
        M_j = np.array(M_j, order='F')
        # print mj, len(mj)
        # print M_j
    else:
        M_j = []
        M_j.append(M)
        tmp_mj = []
        for i in mj:
            M_j.append(len(i))
            for j in i:
                tmp_mj.append(j)
        hyper_join = np.array(tmp_mj, order='F')
        M_j = np.array(M_j, order='F')
        # print mj, len(mj)
        # print M_j
    # tot_hyp = np.sum(M_j)
    # print tot_hyp
    if csr.shape[0] != rhs.shape[0]:
        print "Error!  Mismatch of the matrix and the right hand side sizes!"
        info = 0,0,0,0
        return 0, info
    if np.sum(bl_ar) != csr.shape[0] :
        print "Error! Incorrect block sizes or number of blocks!"
        info = 0,0,0,0
        return 0, info
    if np.sum(hyper_join[:M_j[1]]) != M:
            print np.sum(hyper_join[:M_j[1]]), M
            print "Error!  Incorrect multilevel joining size!"
            info = 0,0,0,0
            return 0, info
    # print hyper_join
    x, mem, f_time, sol_time, res, error = prec.factor_sol(csr.indptr+1, csr.indices+1,
                                         csr.data, rhs,#,level,
                                         eps, r, block_size,bl_ar,
                                         join,hyper_join,M_j)
    # print f_time, sol_time
    info = (f_time, sol_time), mem, res, error
    return x, info

itr = []
def _simple_print(r):
    itr.append(r)
    print len(itr), np.linalg.norm(r)
def _simple_print_cg(xk):
    r = csr.dot(xk)-rhs
    itr.append(r)
    print len(itr), np.linalg.norm(r)
def _report(xk):
    #global csr, rhs
    #print csr.dot(xk) - rhs
    frame = inspect.currentframe().f_back
    itr.append((frame.f_locals['resid'], time()-t1))
    print len(itr), frame.f_locals['resid']


from time import time 
def _count_it(r):
    #xk = r
    #print xk
    frame = inspect.currentframe().f_back
    itr.append((frame.f_locals['resid'], time()-t1))
def _prec_fun(v):
    return prec.solve_ll(v)

def prec_gmres(csr, rhs, eps_it=1e-2, block_size=8, bl_ar=[0], level=6, r=3, restart=15,
               maxiter=100, tol=1e-10, use_prec=True, verbose=1, join=2, mj=[[0]]):
    global itr, t1
    if bl_ar[0] == 0:
        M = csr.shape[0]/block_size
        bl_ar = np.ones(M, order='F')*8
    else:
        M = bl_ar.shape[0]
    N = csr.shape[0]
    #M = N/block_size
    if mj[0][0] == 0:
        M_j = []
        M_j.append(M)
        mj = []
        for i in xrange(level-1):
            if M_j[-1]%join == 0:
                M_j.append(M_j[-1]/join)
                for i in xrange(M_j[-1]):
                    mj.append(join)
            else:
                M_j.append(M_j[-1]/join+1)
                for i in xrange(M_j[-1] - 1):
                    mj.append(join)
                mj.append(M_j[-2]%join)
        hyper_join = np.array(mj, order='F')
        M_j = np.array(M_j, order='F')
        # print mj, len(mj)
        # print M_j
    else:
        M_j = []
        M_j.append(M)
        tmp_mj = []
        for i in mj:
            M_j.append(len(i))
            for j in i:
                tmp_mj.append(j)
        hyper_join = np.array(tmp_mj, order='F')
        M_j = np.array(M_j, order='F')
        # print mj, len(mj)
        # print M_j
    if csr.shape[0] != rhs.shape[0]:
        print "Error!  Mismatch of the matrix and the right hand side sizes!"
        info = 0,0,0,0
        return 0, info
    if np.sum(bl_ar) != csr.shape[0] :
        print "Error! Incorrect block sizes or number of blocks!"
        info = 0,0,0,0
        return 0, info
    if np.sum(hyper_join[:M_j[1]]) != M:
        print np.sum(hyper_join[:M_j[1]]), M
        print "Error!  Incorrect multilevel joining size! Error in level"
        info = 0,0,0,0
        return 0, info
    x0 = np.zeros(N)*1.
    prec_lo = LinearOperator((N, N), _prec_fun, dtype=float)
    prec.free_mem()
    t0 = time()
        #if M_j.shape[0] > 4:
        #print "WHF?!"
        #info = 0,0,0,0
        #return 0, info

    mem, f_time = prec.factor_ll(csr.indptr+1, csr.indices+1, csr.data,
                         eps_it, r, block_size, bl_ar,
                         join,hyper_join,M_j)
    # print "Fact time", f_time
    t_fact = f_time# time() - t0
    t1 = time()
    if verbose == 1:
        if use_prec is True:
            # x1, errr = gmres(csr, rhs, x0, tol, restart, maxiter,
            #                  xtype=0, M=prec_lo, callback=_simple_print)
            x1, errr = cg(csr, rhs, x0, tol, maxiter,
                            xtype=0, M=prec_lo, callback= _report)
        else:
            x1, errr = cg(csr, rhs, x0, tol, maxiter,
                             xtype=0, M=None, callback= _report)
        if (errr == 0):
            t_it = time() - t1
        else:
            t_it = time() - t1
            print "H2-prec r=" + str(r) + " fails to converge!"
            info = t_fact, t_it, mem, itr
            itr = []
            return 0, info
    else:
        if use_prec is True:
            x1, errr = cg(csr, rhs, x0, tol, maxiter,
                             xtype=0, M=prec_lo, callback=_count_it)
        else:
            x1, errr = cg(csr, rhs, x0, tol, maxiter,
                             xtype=0, M=None, callback=_count_it)
        if (errr == 0):
            t_it = time() - t1
        else:
            t_it = time() - t1
            print "H2-prec r=" + str(r) + " fails to converge!"
            info = t_fact, t_it, mem, itr
            itr = []
            return 0, info
    info = t_fact, t_it, mem, itr
    itr = []
    return x1, info
