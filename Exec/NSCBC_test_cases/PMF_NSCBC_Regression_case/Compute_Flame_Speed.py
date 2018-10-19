import os
import subprocess
import numpy as np
import sys
import glob
from itertools import islice
import re

dir = sys.argv[1]
if dir == 'x':
    idir = 1
elif dir == 'y':
    idir = 2
elif dir == 'z':
    idir = 3
else:
    sys.exit("Script requires a single argument that must be x, y or z, to choose slice direction")

numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

listing = sorted(glob.glob('plt*'),key=numericalSort)

TimeSl = 'Flame_Speed.txt'

text_file = open(TimeSl, "w")

for filename in listing:

    print filename

    exeName = './fextract.Linux.gfortran.exe '

    cmdArgs = '-p ' + filename + ' -s ' + 'tempAscii.txt' +' -d ' + str(idir) + ' -v \'density rho_omega_CH4 Y(CH4) x_velocity Temp \''

    callArgs = exeName+cmdArgs

    subprocess.call(callArgs, shell=True)
    
    with open('tempAscii.txt') as fin:
        line = list(islice(fin, 1, 2))
        time = line[0].split()

    x,rho,omega,Y_fuel,xVel,Temp = np.loadtxt('tempAscii.txt', unpack=True, comments='#')

    Sl = -np.trapz(omega,x)/(rho[0]*Y_fuel[0])

    text_file.write('%0.14e \t %f \t %f \n'%(float(time[3]),Sl,-np.trapz(omega,x)))
    
