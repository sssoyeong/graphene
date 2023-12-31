from NanoCore import *

# modeling tool
atoms = io.read_xsf('model.xsf')

# simulation object
sim = s2.Siesta(atoms)

# set simulation options
sim.set_option('kgrid', [1,7,1])         # set #k-point
sim.set_option('kshift', [0.0,0.0,0.0])  # set k-center from gamma
sim.set_option('MixingWt', 0.05)         # adjust mixing weight (density)
sim.set_option('PDOS', 1)                # enable PDOS block
sim.set_option('PDOSE', (-5,5,0.1,1001)) # emin, emax, dE, NE

# run siesta
sim.run(mpi=1, nproc=2)

# get PDOS --> LDOS
z_coords, Z, E = s2.get_pldos(sim, -5, 5, broad=0.05, npoints=1001, label='siesta')

# convert to log10 values + minimum correction to avoid -INF
correct = 10**-5
Z = np.log10(Z+correct)

# generate meshgrid
X, Y = np.meshgrid(np.array(z_coords), E)

# customized figure 
import matplotlib.pyplot as plt
fig1 = plt.figure(figsize=(10,12))
fig11 = fig1.add_subplot(212)
levels = np.linspace(1.01*Z.min(), 0.99*Z.max(), 100)

import pylab as plb
cset = plb.contourf(X,Y,Z, levels)
plb.colorbar(cset)
fig11.axis([0, 20, -5, 5])

# atoms 
fig12 = fig1.add_subplot(311)
Y, Z = atoms.get_2dscatter(plane='yz')
fig12.scatter(Z, Y)
fig12.axis([0, 20, 0, 6])

plt.show()
fig1.savefig('1.png')
