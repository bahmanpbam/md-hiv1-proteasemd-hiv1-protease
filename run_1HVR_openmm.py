from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Step 1: Load and fix the protein (HIV protease)
# ✅ Load cleaned PDB structure
fixer = PDBFixer(filename='1HVR_cleaned.pdb')
fixer.findMissingResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(pH=7.0)


# Step 2: Write fixed PDB
with open('1HVR_fixed.pdb', 'w') as f:
    PDBFile.writeFile(fixer.topology, fixer.positions, f)

# Step 3: Load structure and solvate
pdb = PDBFile('1HVR_fixed.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addSolvent(forcefield, model='tip3p', padding=1.0*nanometer)

# Step 4: Build system and simulation
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,
                                 nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Step 5: Minimize and run simulation
print("Minimizing...")
simulation.minimizeEnergy()

simulation.reporters.append(PDBReporter('trajectory.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, temperature=True,
                                              potentialEnergy=True, step=True))

print("Running MD...")
simulation.step(10000)  # ~20 ps
print("Done.")