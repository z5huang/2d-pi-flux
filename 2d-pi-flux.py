# Time-stamp: <2018-02-15 14:28:13 zshuang>
# 2D pi flux model, Zhang et al 2011, PRB 84 075128

import sys
import os
if'../lib/' not in sys.path:
    sys.path.append('../lib/')
import numpy as np
import math
from scipy import sparse
from numpy import linalg as la
from scipy.sparse import linalg as sla
#from scipy.optimize import minimize as spmin
import scipy.optimize as so
#from npext import *
import npext
import chern_fast
from helper import *
import cmath
import itertools as it
#from scipy.special import binom as binom
import functools
import collections

import pandas as pd
import json

##
# to get all 4-element subsets of an 8-element set, use:
it.combinations(np.arange(8), 4)

def kron_reduce(*ops):
    """ take the kron product of all inputing operators, _in order_ """
    return functools.reduce(np.kron, ops)

class PiFlux:
    """ Single particle level """
    def __init__(self, t=0.5, nx=20, ny=20):
        """ t: next neighbor hopping """
        self.t = t
        self.nx = nx
        self.ny = ny
    def hk(self,kx,ky):
        """ Bloch Hamiltonian """
        h11 = 2*np.cos(kx)
        h12 = (1 + 2*self.t * np.sin(kx)) + (1 - 2*self.t * np.sin(kx)) * np.exp(-1j*ky)
        res = np.array( [[h11,        h12],
                        [h12.conj(), -h11]] )
        return res
    def export_erg(self,fn='erg'):
        par = open(fn+'.par','w')
        print("""
        t = %g
        """%(self.t), file=par)
        par.close()

        dat = open(fn+'.dat', gen_mode(fn + '.dat'))

        labels = 'kx ky E1 E2'
        header = gen_header(labels)
        print(header, file=dat)


        for x,kx in enumerate(np.linspace(0, 2*np.pi, self.nx, endpoint=False)):
            for y,ky in enumerate(np.linspace(0, 2*np.pi, self.ny, endpoint=False)):
                print('\rkx = %g , ky = %g                    '%(kx,ky), end='')
                hh = self.hk(kx,ky)
                eig = la.eigvalsh(hh)
                print('%g\t%g\t%s'%(kx,ky, '\t'.join(['%g'%ee for ee in eig])), file=dat)
            dat.write('\n')
        dat.close()
        
    def memo_u(self):
        """ memoize the wavefunctions on the momentum mesh
        """
        res = np.zeros((self.nx,self.ny,2,2),dtype=np.complex)
        for x,kx in enumerate(np.linspace(0, 2*np.pi, self.nx, endpoint=False)):
            for y,ky in enumerate(np.linspace(0, 2*np.pi, self.ny, endpoint=False)):
                hh = self.hk(kx,ky)
                eig,u = la.eigh(hh)
                res[x,y] = u
        self.u = res
        return res
    def chern(self):
        """ compute the Chern numbers of the Hofstadter model, cf. def hk() """
        uu = self.memo_u()
        return chern_fast.chern(uu, boundary=False)

class PiFlux_gapless:
    # Tight binding, with +1 -1 +1 -1 ... hopping along x, and all +1 hopping along y
    """ Single particle level """
    def __init__(self, nx=20, ny=20):
        """ t: next neighbor hopping """
        self.nx = nx
        self.ny = ny
    def hk(self,kx,ky):
        """ Bloch Hamiltonian """
        h12 = 1 + np.exp(-1j * kx) + np.exp(1j * ky)
        res = np.asarray(  [[0,           h12 ],
                            [h12.conj(),  0   ]] )
        return res
    def export_erg(self,fn='erg_gapless'):
        dat = open(fn+'.dat', gen_mode(fn + '.dat'))

        labels = 'kx ky E1 E2'
        header = gen_header(labels)
        print(header, file=dat)

        for x,kx in enumerate(np.linspace(0, 2*np.pi, self.nx, endpoint=False)):
            for y,ky in enumerate(np.linspace(0, 2*np.pi, self.ny, endpoint=False)):
                print('\rkx = %g , ky = %g                    '%(kx,ky), end='')
                hh = self.hk(kx,ky)
                eig = la.eigvalsh(hh)
                print('%g\t%g\t%s'%(kx,ky, '\t'.join(['%g'%ee for ee in eig])), file=dat)
            dat.write('\n')
        dat.close()
        
    def memo_u(self):
        """ memoize the wavefunctions on the momentum mesh
        """
        res = np.zeros((self.nx,self.ny,2,2),dtype=np.complex)
        for x,kx in enumerate(np.linspace(0, 2*np.pi, self.nx, endpoint=False)):
            for y,ky in enumerate(np.linspace(0, 2*np.pi, self.ny, endpoint=False)):
                hh = self.hk(kx,ky)
                eig,u = la.eigh(hh)
                res[x,y] = u
        self.u = res
        return res
    def chern(self):
        """ compute the Chern numbers of the Hofstadter model, cf. def hk() """
        uu = self.memo_u()
        return chern_fast.chern(uu, boundary=False)

#class PiGutz(PiFlux):
#    """ gutzwiller proj on top of PiFlux """
class Gutz():
    """ Gutzwiller proj """
    @classmethod
    def gen_conf_cls(cls, nsites=8):
        """ generate all half-filling configurations """
        cc = it.combinations(np.arange(nsites), nsites//2)
        res = np.fromiter( it.chain.from_iterable( cc ) , dtype=np.uint8)
        return res.reshape(-1,nsites//2)
    def __init__(self,nx=2,ny=2,ntx=20,nty=20,conf_file=None):
        self.ntx = ntx
        self.nty = nty
        #super().__init__(nx = nx,ny = ny)
        self.nx = nx
        self.ny = ny
        self.gen_conf(conf_file)

    def gen_conf(self,conf_file = None):
        """ generate real space configurations for half filling """

        if conf_file != None:
            self.conf = np.load(conf_file)
        else:
            nsites = self.nx * self.ny * 2
            res = np.fromiter( it.chain.from_iterable( it.combinations(np.arange(nsites), nsites//2) ), dtype=np.uint8).reshape(-1, nsites//2)
            self.conf = res
        
    def u_realspace(self, tx=0, ty=0):
        """ real-space single particle wavefunctions of the *bottom* band

        tx,ty: twisted boundary phases, from 0 to 2pi
        """

        dx,dy = tx/self.nx, ty/self.ny
        #dim = self.nx*self.ny*2
        #res = np.zeros((dim, dim//2),dtype=np.complex)
        res = []

        kx = np.linspace(0, 2*np.pi, self.nx, endpoint=None) + dx
        ky = np.linspace(0, 2*np.pi, self.ny, endpoint=None) + dy
        # matrix of exp(i kx * x), where kx is row index, x is column index
        eikxx = np.exp(1j * np.outer( kx, np.arange(self.nx) ) ) / np.sqrt(self.nx)
        # matrix of exp(i ky * y)
        eikyy = np.exp(1j * np.outer( ky, np.arange(self.ny) ) ) / np.sqrt(self.ny)

        for exx, kxx in zip(eikxx,kx):
            for eyy, kyy in zip(eikyy,ky):
                eig,u = la.eigh(self.hk(kxx,kyy))
                u = u[:,0]
                #u = la.eigh(self.hk(kxx,kyy))[1][:,0]
                u_full = kron_reduce(exx,eyy,u)
                res.append(u_full)
        return np.asarray(res).T
                
    def manybody_wf_1spin(self, tx=0, ty=0):
        """ obtain the many-body wavefunction of a single spin species at half filling, with twisted boundary phases tx and ty """
        uu = self.u_realspace(tx,ty)
        res = np.asarray(
            [ la.det(uu[c]) for c in self.conf ]
            )
        return res

    def export_proj_weight(self, fn='proj-w'):
        """ compute the projected wavefunction of two spin species, onto single occupation sector """

        header = gen_header('tx ty weight')
        with open(fn + '.dat', gen_mode(fn + '.dat')) as dat:
            print(header, file=dat)
            for tx in np.linspace(0, 2*np.pi, self.ntx, endpoint=None):
                for ty in np.linspace(0, 2*np.pi, self.nty, endpoint=None):
                    print('\rtx,ty = %g, %g               '%(tx,ty), end='')
                    s1wf = self.manybody_wf_1spin(tx,ty)
                    # The spin down wavefunction happens to be the spin up
                    # but in reverse spatial configuration order. We are
                    # also ignoring the signature of the combined real
                    # space configuration, which, since it is independent
                    # of tx and ty, is simply a gauge choice.
                    w = la.norm(s1wf * s1wf[::-1])
                    print('%g\t%g\t%g'%(tx,ty,w), file=dat)
                dat.write('\n')
    def memo_manybody_wf_1spin(self):
        """ memoize 1spin wavefunction, i.e., noninteracting """
        dim = self.conf.shape[0]
        res = np.zeros((self.ntx,self.nty,dim), dtype=np.complex)
        for x,tx in enumerate(np.linspace(0, 2*np.pi, self.ntx, endpoint=False)):
            for y,ty in enumerate(np.linspace(0, 2*np.pi, self.nty, endpoint=False)):
                print('\rmemoizing tx, ty = %d, %d           '%(x,y), end='')
                res[x,y] = self.manybody_wf_1spin(tx,ty)
        return res
    def manybody_wf_fermion(self, wf):
        """ convert non-interacting wf to wf with 1 electron per site """
        return wf * wf[:,:,::-1]
    def manybody_wf_boson(self, wf_1spin):
        """ convert non-interacting wf to boson wfs, i.e., 0 or 2 electrons per site """
        return wf * wf
    
    def chern_manybody(self, wf_func=None):
        """ compute the Chern number of a manybody wavefunction, returns the Chern number and the memoized wavefunction """
        if wf_func == None:
            wf_func = self.memo_manybody_wf_2spin
        wf = wf_func()
        return (*chern_fast.chern(wf), wf)

class PiGutz(PiFlux, Gutz):
    """ Pi Flux Gutzwiller """
    def __init__(self, t=0.5, nx=2, ny=2, ntx=20, nty=20, conf_file=None):
        Gutz.__init__(self, nx = nx, ny = ny, ntx = ntx, nty = nty, conf_file = conf_file)
        PiFlux.__init__(self, t=t, nx=nx, ny=ny)

class PiGutz_gapless(PiFlux_gapless, Gutz):
    """ Pi Flux Gutzwiller """
    def __init__(self, nx=2, ny=2, ntx=20, nty=20, conf_file=None):
        Gutz.__init__(self, nx = nx, ny = ny, ntx = ntx, nty = nty, conf_file = conf_file)
        PiFlux.__init__(self, nx=nx, ny=ny)

def export_manybody_chern(wf, fn='berry_curvature'):
    """ write berry curvature to fn, and its berry phase (over ky) to fn2 """
    #wf = np.load(wf_file)
    ntx,nty,dim = wf.shape
    c,curv = chern_fast.chern(wf)
    curv = curv / (2*np.pi)

    par = open(fn + '.par', 'w')
    txt = """
ntx = %d
nty = %d
wf_dim = %d
"""%(ntx,nty,dim)
    print(txt, file=par)
    par.close()

    np.savetxt(fn + '.dat', curv)

def export_manybody_berry(wf, fn='berry_phase'):
    """ write its berry phase (over ky) to fn """
    #wf = np.load(wf_file)
    ntx,nty,dim = wf.shape
    c,curv = chern_fast.chern(wf)
    curv = curv / (2*np.pi)

    par = open(fn + '.par', 'w')
    txt = """
ntx = %d
nty = %d
wf_dim = %d
"""%(ntx,nty,dim)
    print(txt, file=par)
    par.close()

    wf0 = wf[0]
    wf0s = np.roll(wf0, shift=-1, axis=0)
    overlap = np.sum(wf0 * wf0s.conj(), axis=-1)
    berry0 = np.angle(np.prod(overlap))
    dberry = np.average(curv, axis = 1)
    berry = np.cumsum(dberry) / ntx
    berry = [berry0, *(berry0 + berry)]

    dat = open(fn + '.dat', 'w')
    print('\n'.join([ '%g'%b for b in berry ]), file=dat)
    dat.close()

def export_manybody_weight(wf, fn='gutz_weight'):
    """ write normalization to fn """
    #wf = np.load(wf_file)
    weight = np.sum(wf * wf.conj(), axis=-1)
    ntx,nty,dim = wf.shape

    with open(fn + '.par', 'w') as par:
        txt = """
ntx = %d
nty = %d
wf_dim = %d
"""%(ntx,nty,dim)
        print(txt, file=par)

    np.savetxt(fn + '.dat', np.real(weight))

def export_sv_flow(wf, fn='svflow'):
    """ export singular value flow 

    for each kx, obtain singular values along ky
    """
    ntx,nty,dim = wf.shape
    with open(fn + '.par', 'w') as par:
        txt = """
ntx = %d
nty = %d
wf_dim = %d
"""%(ntx,nty,dim)
        print(txt, file=par)


    header = gen_header('n kx s(n,kx)')

    # indexed by [kx, sid]
    sv = np.asarray([ np.sort(la.svd(wf_ky.T, full_matrices = False)[1]) for wf_ky in wf ])
    #np.savetxt(fn + '.dat', sv)
    with open(fn + '.dat', 'w') as dat:
        print(header, file=dat)
        for sid, sval_kx in enumerate(sv.T):
            txt = '\n'.join([ '%d\t%d\t%g'%(sid, kx, sval) for kx,sval in enumerate(sval_kx) ])
            print(txt, file=dat)
            dat.write('\n')
