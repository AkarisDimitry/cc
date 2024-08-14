import os
import shutil
import subprocess
from ase.io import read
from ase.calculators.vasp import Vasp

def auto_correct_path(path_string):
    if not path_string.startswith('/'):
        path_string = '/' + path_string
    if not path_string.endswith('/'):
        path_string = path_string + '/'
    return path_string

def directory_maker(path_input):
    if isinstance(path_input, str):
        path = auto_correct_path(path_input)
        os.makedirs(path, exist_ok=True)
        print(f'Making {path}')
    elif isinstance(path_input, list):
        for i in path_input:
            path = auto_correct_path(i)
            os.makedirs(path, exist_ok=True)
            print(f'Making {path}')

def submit_vasp_calculation(ase_atoms_object, path_to_files, slurm=False, tasks=34, time_cap='01:00:00', partition='genoa'):
    origin = os.getcwd()
    templatefile = os.getenv("VASP_SLURM_TEMPLATE")
    with open(templatefile, "r") as file:
        template = file.read()
    template = template.format(jobname=ase_atoms_object.get_chemical_formula(mode="metal"), ntasks=tasks, partition=partition, time_cap=time_cap)
    submitfile = os.path.join(path_to_files, "submit.sh")
    with open(submitfile, "w") as file:
        file.write(template)
    if slurm:
        os.chdir(path_to_files)
        subprocess.run(['sbatch', 'submit.sh'])
    os.chdir(origin)

def relaunch_vasp_calculation(path):
    if not os.path.exists(os.path.join(path, 'CONTCAR')):
        print(f"No CONTCAR found in {path}. Skipping.")
        return

    atoms = read(os.path.join(path, 'CONTCAR'))
    sub_folder = 'relaunched_vaspsolv'
    work = os.path.join(path, sub_folder)
    directory_maker(work)
    
    # Copy necessary files
    for file in ['WAVECAR', 'CHGCAR']:
        if os.path.exists(os.path.join(path, file)):
            shutil.copy(os.path.join(path, file), work)
    
    os.chdir(work)
    atoms.calc = Vasp(
        xc='RPBE',
        encut=400,
        kspacing=0.03,
        istart=1,  # Start from WAVECAR
        icharg=1,  # Read charge density from CHGCAR
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
        efield=0.1,  # You may want to adjust this value
        ldipol=True,
        idipol=3,
        
        # VASP-sol parameters for implicit solvation in water
        lsol=True,
        eb_k=78.4,
        tau=0,
        lambda_d_k=3.0,
        nc_k=0.0025,
        sigma_k=0.6,
        
        lwave=True,
        lcharg=True,
        npar=4,
        dipol=[0.5, 0.5, 0.4]
    )
    atoms.calc.write_input(atoms=atoms)
    submit_vasp_calculation(ase_atoms_object=atoms, path_to_files=work, slurm=True, tasks=64, time_cap='01:00:00')

def recursive_relaunch(root_path):
    for root, dirs, files in os.walk(root_path):
        if 'CONTCAR' in files:
            print(f"Relaunching calculation in {root}")
            relaunch_vasp_calculation(root)

# Main execution
root_path = '/path/to/your/main/directory'  # Replace with your actual path
recursive_relaunch(root_path)