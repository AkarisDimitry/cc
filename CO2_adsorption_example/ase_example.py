import numpy as np
from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.io import read, write
import os

def run_vasp_calculation(atoms, potential):
    calc = Vasp(
        xc='RPBE',
        encut=400,
        kspacing=0.03,
        istart=0,
        icharg=2,
        prec='Accurate',
        ediff=1e-5,
        nelm=100,
        nelmin=6,
        nsw=0,
        ibrion=-1,
        isif=2,
        ediffg=-0.01,
        ivdw=11,
        
        # Electric field parameters
        efield=potential,  # Applied potential
        ldipol=True,
        idipol=3,
        
        # Constant potential parameters
        bpotim=0.1,
        smass=0,
        nelect=100,  # Adjust this to your system's total number of electrons
        
        # VASP-sol parameters for implicit solvation in water
        lsol=True,
        eb_k=78.4,
        tau=0,
        lambda_d_k=3.0,
        nc_k=0.0025,
        sigma_k=0.6,
        
        lwave=False,
        lcharg=True,
        npar=4,
    )
    
    atoms.calc = calc
    
    energy = atoms.get_potential_energy()
    forces = atoms.get_forces()
    
    return energy, forces

# Read your initial structure
initial_atoms = read('your_structure.xyz')

# Define the range of potentials to apply (in eV/Angstrom)
potentials = np.linspace(-0.5, 0.5, 11)  # 11 points from -0.5 to 0.5

results = []

for i, potential in enumerate(potentials):
    print(f"Running calculation for potential: {potential:.3f} eV/Angstrom")
    
    # Create a directory for this calculation
    dir_name = f"potential_{potential:.3f}"
    os.makedirs(dir_name, exist_ok=True)
    os.chdir(dir_name)
    
    # Run the calculation
    energy, forces = run_vasp_calculation(initial_atoms.copy(), potential)
    
    # Store the results
    results.append({
        'potential': potential,
        'energy': energy,
        'forces': forces
    })
    
    # Write the relaxed structure
    write('relaxed_structure.xyz', initial_atoms)
    
    # Go back to the parent directory
    os.chdir('..')

# Write out the results
with open('potential_sweep_results.txt', 'w') as f:
    for result in results:
        f.write(f"Potential: {result['potential']:.3f} eV/Angstrom\n")
        f.write(f"Energy: {result['energy']:.6f} eV\n")
        f.write("Forces:\n")
        for force in result['forces']:
            f.write(f"{force[0]:.6f} {force[1]:.6f} {force[2]:.6f}\n")
        f.write("\n")

print("Calculations complete. Results written to 'potential_sweep_results.txt'")