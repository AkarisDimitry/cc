from ase import Atoms
from ase.calculators.vasp import Vasp
from ase.io import read
import os
import subprocess
import glob
import numpy as np
import math

def auto_correct_path(path_string):
    # Check if the path string starts with a "/"
    if not path_string.startswith('/'):
        path_string = '/' + path_string
    
    # Check if the path string ends with a "/"
    if not path_string.endswith('/'):
        path_string = path_string + '/'
    return path_string

def directory_maker(path_input):
    input_type = type(path_input)
    if input_type == str:
        path = auto_correct_path(path_input)
        if not os.path.exists(path):
            os.mkdir(path)
            print(f'Making {path}')
        else:
            print(f'{path} exists')
    elif input_type == list:
        for i in path_input:
            path = auto_correct_path(i)
            if not os.path.exists(path):
                os.mkdir(path)
                print(f'Making {path}')
            else:
                print(f'{path} exists')

def submit_vasp_calculation(ase_atoms_object,path_to_files,slurm=False,tasks=32,time_cap='01:00:00', partition = 'genoa'):
    origin = os.getcwd()
    # get template for the submission script
    templatefile = os.getenv("LORENTZ_SLURM_TEMPLATE")
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
    
path = glob.glob('/path/to/CO2_adsorption_example/structures/*/*/*.POSCAR')

for f in path:
    atoms = read(f,format='vasp')
    a,b,c,alpha,beta,gamma = atoms.cell.cellpar()
    ka = round(40/a)
    kb = round(40/b)
    kc = 1
    direct = '/'.join(f.split('/')[0:-1])
    file = f.split('/')[-1].split('.')[0]
    for potential in np.linspace(-0.5, 0.5, 11):
        work = f'{direct}/{file}_potential_{round(potential,2)}eVA/'
        directory_maker(work)
        os.chdir(work)
        atoms.calc = Vasp(xc='RPBE',
                          encut=400,
                          kpts = [ka,kb,kc],
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
                          efield=potential,  # Applied electric field
                          ldipol=True,
                          idipol=3,
                          
                          # Constant potential parameters
                          # bpotim=0.1,
                          # smass=0,
                          # nelect=100,  # Adjust this to your system's total number of electrons
                          
                          # VASP-sol parameters for implicit solvation in water
                          lsol=True,
                          eb_k=78.4,
                          tau=0,
                          lambda_d_k=3.0,
                          nc_k=0.0025,
                          sigma_k=0.6,
                          
                          lwave=True,
                          lcharg=True,
                          npar=4
                         )
        atoms.calc.write_input(atoms = atoms)
        submit_vasp_calculation(ase_atoms_object=atoms, path_to_files=work,slurm=True,tasks=64,time_cap='01:00:00')


