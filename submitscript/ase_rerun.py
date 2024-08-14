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

def submit_vasp_calculation(ase_atoms_object,path_to_files,slurm=False,tasks=34,time_cap='01:00:00', partition = 'genoa'):
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
    
path = glob.glob('/home/ifilot/vasp/workshop/structures/*x*/Ag_CO2*/*eVA/CONTCAR')

for f in path:
    atoms = read(f,format='vasp')
    direct = '/'.join(f.split('/')[0:-1])+'/'
    sub_folder = 'relaunched_vaspsolv'
    work = os.path.join(direct, sub_folder)
    directory_maker(work)
    files = glob.glob(f'{direct}*')
    filename = f.split('/')[-1].split('.')[0]
    for potential in [0.01,0.05,0.1]:
        work = f'{direct}/{filename}_potential_{round(potential,4)}eVA/'
        directory_maker(work)
        os.chdir(work)
        atoms.calc = Vasp(xc='RPBE',
                          encut=400,
                          kspacing =0.3,
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
                          # VASP-sol parameters for implicit solvation in water
                          lsol=True,
                          eb_k=78.4,
                          tau=0,
                          lambda_d_k=3.0,
                          nc_k=0.0025,
                          sigma_k=0.6,
                          lwave=True,
                          npar=4,
                          dipol=[0.5,0.5,0.4]
                         )
        atoms.calc.write_input(atoms = atoms)
        submit_vasp_calculation(ase_atoms_object=atoms, path_to_files=work,slurm=True,tasks=64,time_cap='01:00:00')

for raiz, dirs, archivos in os.walk(directorio_raiz):
    if 'CONTCAR' in archivos and 'WAVECAR' in archivos and 'CHGCAR' in archivos:
        print(f"Procesando directorio: {raiz}")
            
        # Crear nueva carpeta
        nueva_carpeta = os.path.join(raiz, 'relanzado')
        os.makedirs(nueva_carpeta, exist_ok=True)
            
        # Copiar archivos necesarios
        for archivo in ['WAVECAR', 'CHGCAR']:
            shutil.copy(os.path.join(raiz, archivo), nueva_carpeta)
            
        # Leer estructura de CONTCAR
        atoms = read(os.path.join(raiz, 'CONTCAR'))
            
        # Cambiar al nuevo directorio
        os.chdir(nueva_carpeta)

    direct = '/'.join(f.split('/')[0:-1])
    filename = f.split('/')[-1].split('.')[0]
    for potential in [0.01,0.05,0.1]:
        work = f'{direct}/{filename}_potential_{round(potential,4)}eVA/'
        directory_maker(work)
        os.chdir(work)
        atoms.calc = Vasp(xc='RPBE',
                          encut=400,
                          kspacing =0.3,
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
                          
                          lwave=True,
                          npar=4,
                          dipol=[0.5,0.5,0.4]
                         )
        atoms.calc.write_input(atoms = atoms)
        submit_vasp_calculation(ase_atoms_object=atoms, path_to_files=work,slurm=True,tasks=64,time_cap='01:00:00')

