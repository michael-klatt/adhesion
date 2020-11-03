#!/usr/bin/env python

import numpy as np
import math

class sphere:
    """Spherical model of a bacterium with tethers"""

    def __init__(self, R, T, dT):
        # sphere radius
        self.R = R
        # mean tether length
        self.T = T
        # std dev of tether length
        self.dT = dT
        # tethers: xyz coordinates, adhesion, length
        self.tethers = np.empty((0,5))
        # xyz coordinates of patch centers on sphere
        self.centers = np.empty((0,3))
        # mean value of how many uniform tethers
        self.N_pts_mean = 0
        # how many tethers
        self.N_pts = 0

    def uniform_tethers(self, N_pts):
        """Add tethers sampled uniformly on sphere"""
        """        """
        """Input:  scalar N_pts mean number of points"""
        self.N_pts_mean = N_pts
        self.N_pts = np.random.poisson(N_pts)

        t = np.random.uniform(0,2*math.pi,self.N_pts)
        u = np.random.uniform(-1,1,self.N_pts)
        m = np.sqrt(1-u**2)

        x = m*np.cos(t)
        y = m*np.sin(t)
        z = u

        # actual tether length 
        dR = np.random.normal(self.T, self.dT, (self.N_pts,1))

        pts = np.concatenate([x.reshape(-1,1),
                              y.reshape(-1,1),
                              z.reshape(-1,1),
                              np.ones((self.N_pts,1)),
                              dR],axis=1)

        self.tethers = np.vstack([self.tethers,pts])

    def rsa_patch_centers(self, r_rsa, N_trial):
        """Add patch centers using RSA"""
        """        """
        """Input:  scalar r_rsa radius of caps of RSA process"""
        """        scalar N_trial number of cap insertion trials"""
        t = np.random.uniform(0,2*math.pi,N_trial)
        u = np.random.uniform(-1,1,N_trial)
        m = np.sqrt(1-u**2)

        x = m*np.cos(t)
        y = m*np.sin(t)
        z = u

        f = np.concatenate([x.reshape(-1,1),
                            y.reshape(-1,1),
                            z.reshape(-1,1)],axis=1)

        check = pairwise_distances(f) > r_rsa
        check[np.tril_indices(N_trial,k=0)] = True
        keep = np.all(check,axis = 0)

        self.centers = np.vstack([self.centers,f[keep,:]])

    def color_patches(self, adhesion, r_patch):
        """For all patches: change adhesion of tethers in patches"""
        """        """
        """Input:  scalar adhesion new value of adhesion in patch"""
        """        scalar r_patch radius of patch"""
        diff = np.abs(self.tethers[:,:3,None]-self.centers[:,:,None].T)
        D = np.sqrt( (diff**2).sum(1) )

        ids = np.where(D < r_patch)[0]

        self.tethers[ids,3] = adhesion

    def cluster_patches(self, ratio, r_patch):
        """For all patches: change adhesion of tethers in patches"""
        """        """
        """Input:  scalar ratio (> 1) of mean numbers of tethers: patch to background"""
        """        scalar r_patch radius of patch"""
        assert(ratio > 1)

        # How many additional tethers?
        N_t = np.random.poisson((ratio-1) * self.N_pts_mean)

        t = np.random.uniform(0,2*math.pi,N_t)
        u = np.random.uniform(-1,1,N_t)
        m = np.sqrt(1-u**2)

        x = m*np.cos(t)
        y = m*np.sin(t)
        z = u

        # actual tether length 
        dR = np.random.normal(self.T, self.dT, (N_t,1))

        pts = np.concatenate([x.reshape(-1,1),
                              y.reshape(-1,1),
                              z.reshape(-1,1),
                              np.ones((N_t,1)),
                              dR],axis=1)

        # Find tethers in patches and only keep those
        diff = np.abs(pts[:,:3,None]-self.centers[:,:,None].T)
        D = np.sqrt( (diff**2).sum(1) )

        ids = np.where(D < r_patch)[0]

        self.tethers = np.vstack([self.tethers,pts[ids,:]])
        self.N_pts += len(ids)

def pairwise_distances(x):
    """PWD of points (rows of x)"""
    """        """
    """Input:  (N,d) matrix of d-dim coordinates of N points"""
    """Output: (N,N) matrix of pwd"""
    # Trick from stackoverflow.com/questions/28687321
    diff = np.abs(x[:, :, None] - x[:, :, None].T)
    return np.sqrt( (diff**2).sum(1) )

