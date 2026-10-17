"""Independent variable-size lattice and matched production-lag controls.

Existing 24x24 production kernels are left unchanged so their frozen caches
remain valid. These helpers support the supplementary reviewer experiments.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import bmat,coo_matrix,diags,eye

from geometry.lattice import NEIGHBOR_OFFSETS
from model.equations import drift,hill_D,hill_ASC
from model.reduced import drift_DNA
from model.jacobian import jacobian_block


def adjacency(shape,periodic=True):
    rows,cols=shape
    ri,ci=[],[]
    for r in range(rows):
        for c in range(cols):
            for dr,dc in NEIGHBOR_OFFSETS:
                rr,cc=r+dr,c+dc
                if periodic:
                    rr,cc=rr%rows,cc%cols
                if 0<=rr<rows and 0<=cc<cols:
                    ri.append(r*cols+c);ci.append(rr*cols+cc)
    return coo_matrix((np.ones(len(ri)),(ri,ci)),shape=(rows*cols,rows*cols)).tocsr()


def lattice_rhs(t,y,p,W,model="full",notch_multiplier=None):
    ns=4 if model=="full" else 3
    fields=y.reshape(ns,-1)
    D,N=fields[:2]
    f=drift(*fields,p,W@D,W@N) if model=="full" else drift_DNA(*fields,p,W@D,W@N)
    out=np.array(f)
    if notch_multiplier is not None:
        out[1]+=p.lam_N*(notch_multiplier-1)
    return out.ravel()


def lattice_jacobian(y,p,W,model="full"):
    """Exact sparse derivative at arbitrary patterned fields, no symmetry assumption."""
    ns=4 if model=="full" else 3
    D,N,A,*rest=y.reshape(ns,-1)
    sD,sN=W@D,W@N
    if model=="full":S=rest[0]
    else:S=p.lam_ASC/p.d_ASC*hill_ASC(A,p)
    hd,ha=hill_D(S,p),hill_ASC(A,p)
    dhd=p.n_d*hd*(1-hd)/S
    dha=-p.n_ASC*ha*(1-ha)/A
    DD=diags(-p.d_D-p.f_D*sN-p.k_cis*N)
    DN=-diags(p.k_cis*D)-diags(p.f_D*D)@W
    ND=-diags(p.k_cis*N)-diags(p.f_N*N)@W
    NN=diags(-p.d_N-p.f_N*sD-p.k_cis*D)
    AD=diags(p.f_N*N)@W
    AN=diags(p.f_N*sD)
    AA=diags(np.full(len(D),-p.d_A))
    if model=="DNA":
        DA=diags(p.lam_D*dhd*p.lam_ASC/p.d_ASC*dha)
        return bmat([[DD,DN,DA],[ND,NN,None],[AD,AN,AA]],format="csc")
    DS=diags(p.lam_D*dhd)
    SA=diags(p.lam_ASC*dha)
    SS=diags(np.full(len(D),-p.d_ASC))
    return bmat([[DD,DN,None,DS],[ND,NN,None,None],[AD,AN,AA,None],[None,None,SA,SS]],format="csc")


def production_lag_block(mu,fp,p):
    """D,N,A,U; U relaxes toward the matched composite Delta production."""
    J=jacobian_block(mu,fp,p)
    kappa=J[0,3]
    out=np.zeros((4,4))
    out[:3,:3]=J[:3,:3]
    out[0,3]=1
    out[3,2]=kappa*J[3,2]
    out[3,3]=-p.d_ASC
    return out


def lag_similarity(fp,p):
    return np.diag([1.,1.,1.,jacobian_block(0,fp,p)[0,3]])


def random_initial(fp,p,shape,seed,amp=1e-4,free_S=False):
    rng=np.random.Generator(np.random.PCG64(seed))
    xi=rng.uniform(-1,1,(4,*shape))
    xi-=xi.mean(axis=(1,2),keepdims=True)
    xi/=np.max(np.abs(xi),axis=(1,2),keepdims=True)
    fields=fp[:,None,None]*(1+amp*xi)
    if not free_S:
        fields[3]=p.lam_ASC/p.d_ASC*hill_ASC(fields[2],p)
    return fields
