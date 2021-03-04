# adhesion
Simulate patchy spheres and estimate the adhesion forces at sinusoidal surfaces.

* Author: Michael A. Klatt
* E-Mail: software@mklatt.org
* Website: https://mklatt.org
* License: GNU GPLv3

## 

The code was written and used to simulate a geometric model of bacteria with adhesive molecules in the cell envelope in the following paper:

Christian Spengler, Bernhard A. Glatz, Erik Maikranz, Markus Bischoff, Michael Andreas Klatt, Ludger Santen, Andreas Fery, Karin Jacobs. "The adhesion capability of S. aureus cells is heterogeneously distributed over the cell envelope" bioRxiv 2021.01.05.425282
https://doi.org/10.1101/2021.01.05.425282 

If you use the code in your research, please cite the paper.

## ToDo

* [ ] Revise flexible tethers
* [ ] Speed up: numba function: wave.parallel_interpolated

## Help output

```
usage: adhesion.py [-h] [-p PREFIX] [--n_runs N_RUNS] [--seed SEED]
                   [-A AMPLITUDE] [-l WAVELENGTH] [--dx DX] [-R RADIUS]
                   [-T TETHER] [-dT DT] [-P PATCH] [--r_rsa R_RSA]
                   [--n_rsa N_RSA] [-n N] [--adhesion ADHESION]
                   [--cluster CLUSTER] [--rigid]

Simulate patchy spheres and estimate their adhesion (all lengths in mu)

optional arguments:
  -h, --help            show this help message and exit

Simulation:
  -p PREFIX, --prefix PREFIX
                        Prefix of output files (default: ../patchy-data/)
  --n_runs N_RUNS       Number of simulation runs (default: 14)
  --seed SEED           Seed of simulation runs (0: no seed) (default: 0)

Profile:
  -A AMPLITUDE, --amplitude AMPLITUDE
                        Amplitude of surface profile (default: 0.3)
  -l WAVELENGTH, --wavelength WAVELENGTH
                        Wavelength of surface profile (default: 2.75)
  --dx DX               Step-size in x (default: 0.01)
  --rigid               Are the tethers rigid? (default: True)

Sphere:
  -R RADIUS, --radius RADIUS
                        Radius of sphere (default: 0.5)
  -T TETHER, --tether TETHER
                        Mean tether length (default: 0.05)
  -dT DT                Std. dev. of tether length (default: 0.001)

Patches:
  -P PATCH, --patch PATCH
                        Radius of patch (default: 0.125)
  --r_rsa R_RSA         Radius of RSA cap (default: 0.85)
  --n_rsa N_RSA         Mean number of RSA cap insertion trials (default:
                        1000)

Tethers:
  -n N                  Mean num. of unif. distr. tethers (default: 10000)
  --adhesion ADHESION   Increased adhesion in patches (default: 1.0)
  --cluster CLUSTER     Ratio of mean numbers of tethers: patch to background
                        (default: 0.0)
```

