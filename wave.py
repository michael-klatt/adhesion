#!/usr/bin/env python

import numpy as np
import math

from scipy import optimize

import matplotlib.pyplot as plt

class profile:
    """Wave surface profile"""

    def __init__(self, A, l, dx):
        # scalar amplitude A
        self.A = A
        # scalar wavelength l
        self.l = l
        self.hl = l*0.5
        # scalar angular frequency o (omega)
        self.o  = 2*math.pi/l
        # step-size in x
        self.dx = dx
        # x coordinates of profile (N vector)
        self.x = np.arange(-self.hl,self.hl,self.dx)
        self.N = len(self.x)
        # z coordinates of profile (N vector)
        self.z = self.wave(self.x)
        # derivative of profile (N vector)
        self.D = self.D_wave(self.x)

        # Auxiliary variables for interpolated parallel curve
        # x coordinates of profile (N vector)
        self.px = np.arange(-self.l*0.75,self.l*0.75,self.dx)
        self.pN = len(self.px)
        # z coordinates of profile (N vector)
        self.pz = self.wave(self.px)
        # derivative of profile (N vector)
        self.pD = self.D_wave(self.px)
        # distance variables
        self.tmp_a = np.sqrt(1+self.pD**2)
        self.tmp_b = self.pD/self.tmp_a


    def wave(self, x):
        """Surface profile"""
        return self.A*np.cos(self.o * x)

    def D_wave(self, x):
        """Surface profile"""
        return -self.A*self.o*np.sin(self.o * x)

    def parallel_wave(self, d):
        """Parallel surface: (x,z)-coordinates"""
        """        """
        """Input: scalar d distance of parallel surface """
        """Output: (N,2) matrix of (x,z) coordinates of centers of sphere"""
        # Coordinates of parallel curve
        x_d = self.x - d*self.D/np.sqrt(1+self.D**2)
        z_d = self.z + d       /np.sqrt(1+self.D**2)

        # Potential cusp with non-unique projection
        idx = x_d//self.hl == self.x//self.hl

        out = np.empty((np.sum(idx),2))
        out[:,0] = x_d[idx]
        out[:,1] = z_d[idx]
        return out

    def parallel_interpolated(self, x, d):
        """Parallel surface: (x,z)-coordinates"""
        """        """
        """Input:  N vector of x coordinates"""
        """        N vector distances of parallel surface """
        """Output: N vector of z coordinate of parallel surface"""
        # Coordinates of parallel curve
        x_d = self.px[:,None] - d[None,:] * self.tmp_b[:,None] 
        z_d = self.pz[:,None] + d[None,:] / self.tmp_a[:,None]

        # Potential cusp with non-unique projection
        idx = x_d//self.hl == (self.px//self.hl)[:,None]

        x_p = x_d[idx].reshape(self.pN,-1)
        z_p = z_d[idx].reshape(self.pN,-1)

        bins = x[None,:] < x_p
        digitize = np.argmax(bins,axis=0)
        bigitize = digitize-1
        i = 10
        one = z_p[digitize,np.arange(len(digitize))]
        two = z_p[bigitize,np.arange(len(bigitize))]

        return (one + two)*0.5

    def distance(self, xyz):
        """Distance of points to surface (slow)"""
        """        """
        """Input:  (N,3) matrix of coordinates of points"""
        """Output: N vector of distances"""
        N = np.shape(xyz)[0]
        dist = np.empty(N)
        for i in range(N):
            x = xyz[i,0]
            z = xyz[i,2]

            limits = (x//self.hl)*np.array([0,self.hl])
            root = optimize.root_scalar(root_finding, x0=x,
                                        args=(x,z,self),
                                        method='brentq',
                                        bracket=limits)
            xs = root.root
            zs = self.wave(xs)

            dist[i] = np.sqrt((xs-x)**2+(zs-z)**2)

        return dist

    def show_waves(self,d):
        plt.plot(self.x,self.z,'-')

        centers = self.parallel_wave(d)
        plt.plot(centers[:,0],centers[:,1],'-')

        plt.show()

def root_finding(xs, x, z, surface):
    """Zero of function is x-coordinate of closest point"""
    return x - xs + (z - surface.wave(xs))* surface.D_wave(xs)



