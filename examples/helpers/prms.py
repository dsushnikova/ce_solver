def kdt_prm(msh_coo, B=16):
    import numpy as np
    import scipy as sp
    from sklearn.neighbors import KDTree
    #msh_coo = mesh.coordinates()
    tree = KDTree(msh_coo, leaf_size=2)
    prm = np.array(tree.idx_array)
    N = msh_coo.shape[0]
    nbl = N/B
    bl_ar = np.array([B]*nbl)
    if N%B != 0:
        bl_ar[-1] += N%B
    return prm, bl_ar


def metis_prm(csr, B = 16):
    import metis
    import networkx as nx
    import  pyamg
    import numpy as np
    from time import time
    bl_ar = []
    N = csr.shape[0]
    nparts = N/B
    G = nx.from_scipy_sparse_matrix(csr)
    n, l = metis.part_graph(G, nparts=nparts, recursive=True)
    prm = []
    tmp = []
    for i in xrange(nparts): tmp.append([])
    k = 0
    for i in l:
        tmp[i].append(k)
        k += 1
    k = 0
    for i in tmp:
        bl_ar.append(len(i))
        k += 1
        for j in i:
            prm.append(j)
    return prm, np.array(bl_ar)







