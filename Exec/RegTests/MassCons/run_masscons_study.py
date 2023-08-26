import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

inputs = ['mol-1', 'mol-2', 'plm', 'ppm', 'isothermal', 'isothermal-whydro']

for iname in inputs:
    os.system('rm -r datlog plt*0')
    os.system('./PeleC2d.gnu.MPI.ex masscons-{}.inp max_step=200'.format(iname))
    data = pd.read_fwf('datlog')
    plt.figure('mass')
    plt.plot(data.index, (data['mass']-data['mass'][0])/data['mass'][0], label=iname)

    if iname is not 'isothermal-whydro':
        plt.figure('energy')
        plt.plot(data.index, (data['rho_E']-data['rho_E'][0])/data['rho_E'][0], label = iname)

plt.figure('mass')
plt.xlabel('Timestep')
plt.ylabel('Relative Mass Change')
plt.ylim([-1e-14, 1e-14])
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('figure_conservation_mass.png')
plt.clf()
plt.close()

plt.figure('energy')
plt.xlabel('Timestep')
plt.ylabel('Relative Energy Change')
plt.ylim([-1e-14, 1e-14])
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig('figure_conservation_energy.png')
plt.clf()
plt.close()
