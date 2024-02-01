import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(current_dir)

from NexaBind.app import *
from NexaBind import *
from NexaBind.unit import *


pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
#integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system)
simulation.context
simulation.minimizeEnergy
simulation.reporters
simulation.reporters
simulation.step
