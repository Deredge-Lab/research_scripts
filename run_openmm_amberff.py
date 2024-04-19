from __future__ import print_function

import simtk.openmm.app as app 
import simtk.openmm as mm
import simtk.unit as u
import numpy as np

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

# Add input files (Amber)
prmtop = app.AmberPrmtopFile('CaM_apo.prmtop')
inpcrd = app.AmberInpcrdFile('CaM_apo.inpcrd')

# Set the number of steps for each stage
min_res_max_cycles = 5000                  # Max minimization cycles w/ restraints (default = 5000 cycles) ~CAUTION: Overminimizing can cause the protein to get stuck in an energy minimum
min_nores_max_cycles = 5000                # Max minimization cycles w/o restraints (default = 5000 cycles)
min_tol = 10.0*u.kilojoule_per_mole        # Energy threshold to stop minimization (kJ/mol) (default = 10 kJ/mol)
heating_steps = 100                        # Number of heating steps
heating_interval = 3*u.kelvin              # Increase the temp by 'x' degrees K per 1000 steps
equil1_steps = 100000                      # Equilibrate with restraints for stages 1-5 (default = 0.2 ns per stage) 
equil2_steps = 100000                      
equil3_steps = 100000                      
equil4_steps = 100000                      
equil5_steps = 100000                     
equil6_steps = 500000                      # Equilibrate w/o restraints (default = 1 ns)
prod_steps = 100000000                     # Run production w/o restraints (default = 200 ns)  
dt = 2.0*u.femtoseconds                    # Timestep (default = 2 fs if bonds containing H-atoms are constrained)

# Set restraints for protein heavy atoms
min1_k = 500.0                             # Apply strong restraints to prevent movement of protein BB and SC atoms (default k = 500 kJ/mol/A**2) 
min2_k = 0.0                               # Remove restraints
heating_k = 10.0                           # Apply weak restraints to allow some movement (default k = 10.0 kJ/mol/A**2)
equil1_k = 10.0                            # Begin equilibration, gradually reducing restraints over each stage
equil2_k = 5.0                             # Reduce restraints (default k = 5.0 kJ/mol/A**2)
equil3_k = 2.5                             # Reduce restraints (default k = 2.5 kJ/mol/A**2)
equil4_k = 1.0                             # Reduce restraints (default k = 1.0 kJ/mol/A**2)
equil5_k = 0.1                             # Reduce restraints (default k = 0.1 kJ/mol/A**2)
equil6_k = 0.0                             # Remove restraints 
prod_k = 0.0                    

# Set simulation parameters
fric_coeff = 1.0/u.picoseconds             # Particle collision frequency, gamma (default = 1 per ps)
T_init = 0.0*u.kelvin                      # Initial temperature (default = 0K)
T_fin = 300.0*u.kelvin                     # Final temerature after heating stage (default = 300K)
P = 1.0*u.bar                              # Pressure to be maintained by MC barostat (default = 1 bar)            
update_boxvol_step_interval = 25           # Attempt to update the box volume every 'x' number of steps (default = 25 steps)
add_constraints = app.HBonds               # Constrain the length of hydrogen-containing bonds 
nb_method = app.PME                        # Use the Particle-mesh Ewald algorithm to model nonbonded interactions
nb_cutoff_dist = 0.9*u.nanometers          # Apply a cutoff distance for nonbonded interactions (default = 9 A)
ewald_error_tol = 0.0005                   # Set the error tolerance for PME nonbonded interactions (default = 0.0005)
integrator_constraint_tol = 0.00001        # Set the constraint tolerance for the Langevin integrator (default = 0.00001)
#rigidWater                                # Model waters as solid triangles (rigid) vs. allow them to be flexible (default = True)
#removeCMMotion = False                    # Remove the center of mass motion at every time step to prevent the system from drifting (default = True, unless Langevin dynamics are used)    

# System setup
print('Constructing the system...')
system = prmtop.createSystem(nonbondedMethod=nb_method, nonbondedCutoff=nb_cutoff_dist, constraints=add_constraints, 
    ewaldErrorTolerance=ewald_error_tol, rigidWater=True, removeCMMotion=False)        

print('Adding a harmonic force to restrain protein heavy atoms...s')
if system.usesPeriodicBoundaryConditions():
    energy_expression = '0.5*k*periodicdistance(x, y, z, x0, y0, z0)^2' # periodic distance
    print('System uses PBC') 
else:
    energy_expression = '(K/2)*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)' # non-periodic distance
    print('System does not use PBC') 
restraint_force = mm.CustomExternalForce(energy_expression) 
restraint_force.addGlobalParameter('k', 0*u.kilocalorie_per_mole/u.angstrom**2)
restraint_force.addPerParticleParameter("x0")
restraint_force.addPerParticleParameter("y0")
restraint_force.addPerParticleParameter("z0")
atom_list = list(prmtop.topology.atoms())
for i, atom_crd in enumerate(inpcrd.positions):  
    atoms = atom_list[i] 
    if atoms.name not in ('H', 'HOH', 'K+', 'Cl-'): 
        restraint_force.addParticle(i, atom_crd.value_in_unit(u.nanometers))
system.addForce(restraint_force)

print('Adding a MC barostat with isotropic position scaling...') 
barostat = mm.MonteCarloBarostat(P, T_fin, update_boxvol_step_interval)
barostat.usesPeriodicBoundaryConditions=True
system.addForce(barostat)

print('Constructing a Langevin thermostat')
integrator = mm.LangevinIntegrator(T_init, fric_coeff, dt)   
integrator.setConstraintTolerance(integrator_constraint_tol)   

# Select a GPU to run MD - check that it's available first (nvidia-smi)
platform = mm.Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision': 'mixed', 'CudaDeviceIndex': '0'}

print('Precision format:', properties['CudaPrecision'])
print('GPU:', properties['CudaDeviceIndex'])

print('Constructing the simulation...')
sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)

print('Positioning atoms...') 
sim.context.setPositions(inpcrd.positions)

print('Establishing periodic boundary conditions...')
if inpcrd.boxVectors is not None:
    sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

k0 = sim.context.getParameter('k')
print('Default restraint force:', k0, 'kcal/mol/A**2')
p0 = sim.context.getParameter(barostat.Pressure())
print('Default pressure:', p0, 'bar')

print('Initial system energy: ')
print(sim.context.getState(getEnergy=True).getPotentialEnergy())  

# Minimization Stage 1 (NVT)
print('Minimizing with restraints (stage 1)...')
# openMM uses the L-BFGS minimization algorithm (similar in principle to CG, but with faster convergence)
# tleap tends to leave a large gap near the protein surface and box edges during solvation to avoid violating the laws of physics (e.g., due to spatial overlap between protein and water molecules)
# Briefly equilibrating the system in the NVT ensemble often helps to fill in these gaps so "vacuum bubbles" don't destabilize the system
sim.context.setParameter(barostat.Pressure(), 0.0)
sim.context.setParameter(barostat.Temperature(), 0.0)
p1 = sim.context.getParameter(barostat.Pressure())
p2 = sim.context.getParameter(barostat.Temperature())
print('NVT ensemble (pressure =', p1, 'bar')
sim.context.setParameter('k', min1_k) 
k1 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k1, 'kcal/mol/A**2')
with open('min1.txt', 'w') as f: 
    for i in range(min_res_max_cycles):
        f.write(str(i) + '\t' + str(sim.context.getState(getEnergy=True, enforcePeriodicBox=True).getPotentialEnergy())  + '\n')   
        sim.minimizeEnergy(tolerance=min_tol, maxIterations=min_res_max_cycles)
f.close()
    
# Minimization Stage 2 (NVT) 
print('Minimizing without restraints (stage 2)...')
# In this stage, the restraints are removed so the protein heavy atoms can relax
sim.context.setParameter('k', min2_k) 
k2 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k2, 'kcal/mol/A**2')      
with open('min2.txt', 'w') as f: 
    for i in range(min_nores_max_cycles):
        f.write(str(i) + '\t' + str(sim.context.getState(getEnergy=True, enforcePeriodicBox=True).getPotentialEnergy())  + '\n')   
        sim.minimizeEnergy(tolerance=min_tol, maxIterations=min_res_max_cycles)
f.close()
min2_positions = sim.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()  
app.PDBFile.writeFile(prmtop.topology, min2_positions, open('CaM_apo_min.pdb', 'w'))   

# Generation of particle velocities
# If T_init = 0, velocities are calculated from the forces; if T_init > 0, velocities are randomly chosen from a Boltzmann distr. of values for that temp
print('Setting initial particle velocities...')
sim.context.setVelocitiesToTemperature(T_init)

# Heating (NVT) 
# The canonical (NVT) ensemble is used throughout minimization and heating b/c the system needs to be thermally equilibrated first
# Otherwise, problems can occur due to barostat overcorrections (e.g., bonds containing H atoms aren't properly constrained, water enters hydrophobic pockets of the protein, etc.)
print('Heating...')
sim.context.setParameter('k', heating_k)
k3 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k3, 'kcal/mol/A**2')
with open('heat.txt', 'w') as f: 
    for i in range(heating_steps):      
        integrator.setTemperature(T_fin-(heating_interval*(heating_steps-i)))       
        f.write(str(i) + '\t' + str(integrator.getTemperature()) + '\t' + str(sim.context.getState(getEnergy=True, enforcePeriodicBox=True).getPotentialEnergy())  + '\n')  
        sim.step(1000)
f.close()

# Equilibration Stage 1 (NPT)  
# Now that the temperature of the system is equilibrated, the density of the water molecules needs to be equilibrated by adjusting the box volume
print('Switching to NPT ensemble')
integrator.setTemperature(T_fin)
sim.context.setParameter(barostat.Temperature(), T_fin)
sim.context.setParameter(barostat.Pressure(), P)
print('Equilibrating with restraints (stage 1)...')
sim.context.setParameter('k', equil1_k) 
k4 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k4, 'kcal/mol/A**2')
# Convert from .dcd to .nc using mdconvert, if needed
sim.reporters.append(DCDReporter('equil1_traj.dcd', 1000))
sim.step(equil1_steps)

# Equilibration Stage 2 (NPT)  
print('Equilibrating (stage 2)...')
sim.context.setParameter('k', equil2_k) 
k5 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k5, 'kcal/mol/A**2')
sim.reporters.append(DCDReporter('equil2_traj.dcd', 1000))
sim.step(equil2_steps)
    
# Equilibration Stage 3 (NPT) -
print('Equilibrating (stage 3)...')
sim.context.setParameter('k', equil3_k) 
k6 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k6, 'kcal/mol/A**2')
sim.reporters.append(DCDReporter('equil3_traj.dcd', 1000))
sim.step(equil3_steps)

# Equilibration Stage 4 (NPT) 
print('Equilibrating (stage 4)...')
sim.context.setParameter('k', equil4_k)
k7 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k7, 'kcal/mol/A**2')
sim.reporters.append(DCDReporter('equil4_traj.dcd', 1000))
sim.step(equil4_steps)

# Equilibration Stage 5 (NPT) 
print('Equilibrating (stage 5)...')
sim.context.setParameter('k', equil5_k)
k8 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k8, 'kcal/mol/A**2')
sim.reporters.append(DCDReporter('equil5_traj.dcd', 1000))
sim.step(equil5_steps)

# Equilibration Stage 6 (NPT) 
print('Equilibrating without restraints (stage 6)...')
sim.context.setParameter('k', equil6_k)
k9 = sim.context.getParameter('k')
print('Positional restraint force constant (k):', k9, 'kcal/mol/A**2')
sim.reporters.append(DCDReporter('equil6_traj.dcd', 1000))
sim.step(equil6_steps) 

# Production (NPT) 
print('Starting production...')
sim.context.setParameter('k', prod_k)
sim.reporters.append(app.StateDataReporter('prod.txt', 1000, step=True, time=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True, 
    temperature=True, volume=True, density=True, progress=True, remainingTime=True, speed=True, totalSteps=prod_steps, separator='\t'))  
# The checkpoint and trajectory writeout frequencies should be equal 
sim.reporters.append(DCDReporter('prod_traj.dcd', 10000))
sim.reporters.append(CheckpointReporter('prod.chk', 10000))
sim.step(prod_steps) 
