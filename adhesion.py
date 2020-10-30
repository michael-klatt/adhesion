#!/usr/bin/env python

import numpy as np
import numba

import argparse
import os

from wave   import profile
from patchy import sphere

def parse():
    """Parse parameters from command line"""
    d = "Simulate patchy spheres and estimate their adhesion (all lengths in mu)"
    parser = argparse.ArgumentParser(description=d,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    group_1 = parser.add_argument_group(title='Simulation')
    group_2 = parser.add_argument_group(title='Profile')
    group_3 = parser.add_argument_group(title='Sphere')
    group_4 = parser.add_argument_group(title='Patches')
    group_5 = parser.add_argument_group(title='Tethers')

    # Simulation
    defpref = '../patchy-data/'
    group_1.add_argument('-p', '--prefix',     type=str, action='store',
                         default=defpref,      help='Prefix of output files')
    group_1.add_argument('--n_runs',           type=int, action='store',
                         default=10,           help='Number of simulation runs')
    # Profile
    group_2.add_argument('-A', '--amplitude',  type=float, action='store',
                         default=0.30,         help='Amplitude of surface profile')
    group_2.add_argument('-l', '--wavelength', type=float, action='store',
                         default=2.75,         help='Wavelength of surface profile')
    group_2.add_argument('--dx',               type=float, action='store',
                         default=0.01,         help='Step-size in x')
    # Sphere
    group_3.add_argument('-R', '--radius',     type=float, action='store',
                         default=0.50,         help='Radius of sphere')
    group_3.add_argument('-T', '--tether',     type=float, action='store',
                         default=0.05,         help='Mean tether length')
    group_3.add_argument('-dT',                type=float, action='store',
                         default=0.001,        help='Std. dev. of tether length')
    # Patches
    group_4.add_argument('-P', '--patch',      type=float, action='store',
                         default=0.125,        help='Radius of patch')
    group_4.add_argument('--r_rsa',            type=float, action='store',
                         default=0.75,         help='Radius of RSA cap')
    group_4.add_argument('--n_rsa',            type=int, action='store',
                         default=1000,         help='Number of RSA cap insertion trials')
    # Tethers
    group_5.add_argument('-n',                 type=int, action='store',
                         default=10000,        help='Mean num. of unif. distr. tethers')
    group_5.add_argument('--adhesion',         type=float, action='store',
                         default=5.0,          help='Increased adhesion in patches')
    group_5.add_argument('--n_pp',             type=float, action='store',
                         default=50.0,         help='Mean number of additional tethers per patch')
    group_2.add_argument('--rigid',            action="store_true",
                         default=False,        help='Are the tethers rigid?')
    return parser.parse_args()

def main(args):
    """Main function"""
    os.makedirs(args.prefix, exist_ok=True)

    for run in range(args.n_runs):
        main_run(args, args.prefix, run)

def main_run(args, prefix, run):
    # Patchy Sphere
    sph = sphere(args.radius, args.tether, args.dT)
    #
    N_thethers = np.random.poisson(args.n)
    sph.uniform_tethers(N_thethers)
    #
    # Patches
    sph.rsa_patch_centers(args.r_rsa, args.n_rsa)
    if args.adhesion != 1.0:
        sph.color_patches(args.adhesion, args.patch)
    if args.n_pp > 0:
        sph.cluster_patches(args.n_pp, args.patch)

    # Wave profile
    wav = profile(args.amplitude, args.wavelength, args.dx)
    # wav.show_waves(sph.R)
    centers = wav.parallel_wave(sph.R)

    # Fractions of binding tethers
    if args.rigid:
        frac = np.apply_along_axis(fractions_rigid,
                                   1, centers, sph, wav)
    else:
        frac = np.apply_along_axis(fractions_flexible,
                                   1, centers, sph, wav)

    out = np.column_stack((centers[:,0],frac))

    np.savetxt(prefix + 'sample-' + str(run) + '-fractions.dat', out)
    np.savetxt(prefix + 'sample-' + str(run) + '-points.dat', sph.tethers)
    with open(prefix + 'samples.log','w') as file:
       for arg, value in sorted(vars(args).items()):
        file.write("%s: %r\n"%(arg, value))

# ===================================================================

def fractions_rigid(c, sphere, profile):
    """Fraction of rigid tethers that hit surface"""
    """        """
    """Input:  [x,z] array of sphere center"""
    """        patchy.sphere"""
    """        wave.profile"""
    """Output: fraction of tethers that hit surface"""
    L = (sphere.R + sphere.tethers[:,4])
    tips = L[:,None] * sphere.tethers[:,:3]
    tips[:,0] += c[0] # x coord
    tips[:,2] += c[1] # z coord

    hitting = profile.wave(tips[:,0]) > tips[:,2]

    norm = np.sum(sphere.tethers[:,3])
    frac = np.sum(hitting * sphere.tethers[:,3])/norm

    return frac

def fractions_flexible(c, sphere, profile):
    """Fraction of flexible tethers that hit surface"""
    """        """
    """Input:  [x,z] array of sphere center"""
    """        patchy.sphere"""
    """        wave.profile"""
    """Output: fraction of tethers that reach the surface"""
    L = sphere.tethers[:,4]
    start = sphere.R * sphere.tethers[:,:3]
    start[:,0] += c[0] # x coord
    start[:,2] += c[1] # z coord

    print(c[0])

    # precise but slow distance of tethers to surface
    # hitting = profile.distance(start) < L

    hitting = parallel_interpolated(start[:,0],L,
                  profile.px, profile.pz,
                  profile.tmp_b, profile.tmp_a, 
                  profile.hl, profile.pN) > start[:,2]

    norm = np.sum(sphere.tethers[:,3])
    frac = np.sum(hitting * sphere.tethers[:,3])/norm

    return frac

# Problem: unsupported array index type none in Tuple
# https://stackoverflow.com/a/57916110/13995159
@numba.njit()
def parallel_interpolated(x, d,
                          px, pz, tmp_b, tmp_a, hl, pN):
    """Parallel surface: (x,z)-coordinates"""
    """        """
    """Input:  N vector of x coordinates"""
    """        N vector distances of parallel surface """
    """Output: N vector of z coordinate of parallel surface"""
    # Coordinates of parallel curve
    x_d = px[:,None] - d[None,:] * tmp_b[:,None] 
    z_d = pz[:,None] + d[None,:] / tmp_a[:,None]

    # Potential cusp with non-unique projection
    idx = x_d//hl == (px//hl)[:,None]

    x_p = x_d[idx].reshape(pN,-1)
    z_p = z_d[idx].reshape(pN,-1)

    bins = x[None,:] < x_p
    digitize = np.argmax(bins,axis=0)
    bigitize = digitize-1
    i = 10
    one = z_p[digitize,np.arange(len(digitize))]
    two = z_p[bigitize,np.arange(len(bigitize))]

    return (one + two)*0.5

# ===================================================================

if __name__ == '__main__':
    args = parse()
    main(args)

