from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.io import read

# Read your structure (replace 'your_structure.xyz' with your actual structure file)
atoms = read('your_structure.xyz')

# Create a VASP calculator
calc = Vasp(
    # General VASP parameters
    xc='RPBE',        # Changed to RPBE
    encut=400,        # Confirmed cutoff energy of 400 eV
    istart=0,
    icharg=2,
    prec='Accurate',
    ediff=1e-5,
    nelm=100,
    nelmin=6,
    
    # K-point density
    kspacing=0.03,    # This corresponds to KDENS=0.03
    
    # Ionic relaxation parameters
    nsw=0,
    ibrion=-1,
    isif=2,
    ediffg=-0.01,
    
    # DFT-D3 correction
    ivdw=11,
    
    # Electric field parameters
    efield=0.1,
    ldipol=True,
    idipol=3,
    
    # Constant potential parameters
    bpotim=0.1,
    smass=0,
    nelect=100,  # Adjust this to your system's total number of electrons
    
    # Output options
    lwave=False,
    lcharg=True,
    
    # Parallelization
    npar=4,
    
    # Additional ASE-specific parameters
    directory='vasp_run',
    txt='vasp.out'
)

# Set the calculator for the atoms
atoms.calc = calc

# Run the calculation
atoms.get_potential_energy()

print(atoms)
print(atoms.get_potential_energy())