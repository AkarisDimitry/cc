from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.io import read
import os
import subprocess

def submit_vasp_calculation(ase_atoms_object,path_to_files,slurm=False,tasks=32,time_cap='01:00:00', partition = 'genoa'):
    origin = os.getcwd()
    # get template for the submission script
    templatefile = os.getenv("VASP_SLURM_TEMPLATE")
    with open(templatefile, "r") as file:
        template = file.read()
    # add the chemical formula of the material as jobname in the submission script
    template = template.format(jobname=ase_atoms_object.get_chemical_formula().format("metal"),ntasks=tasks,partition=partition,time_cap=time_cap)
        # write the submission script from template
    submitfile = os.path.join(path_to_files,"submit.sh")
    if not os.path.exists(submitfile):
        with open(submitfile, "w") as file:
            file.write(template)
    # slurm the script (if it exists)
    if os.path.exists(submitfile) and slurm:
        os.chdir(path_to_files)
        subprocess.run(['sbatch','submit.sh'])
    os.chdir(origin)
    
# Read your structure (replace 'your_structure.xyz' with your actual structure file)
atoms = read('your_structure.cif')

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
atoms.calc.write_input(atoms=atoms)
submit_vasp_calculation(ase_atoms_object=atoms,path_to_files=os.get_cwd(),slurm=False,tasks=32,time_cap='01:00:00', partition = 'genoa')
