import os
import shutil
from ase.io import read
from ase.calculators.vasp import Vasp

def run_vasp_calculation(atoms, potential):
    calc = Vasp(
        xc='RPBE',
        encut=400,
        kspacing=0.3,
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
        dipol = [0.5,0.5,0.4],
        
       
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
    atoms.get_potential_energy()  # This will run the calculation

def relanzar_calculos(directorio_raiz):
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
            
            # Relanzar cálculo
            potential = 0.0  # Ajusta este valor según sea necesario
            run_vasp_calculation(atoms, potential)
            
            print(f"Cálculo relanzado en: {nueva_carpeta}")

# Uso del script
directorio_raiz = '/ruta/a/tu/directorio/principal'
relanzar_calculos(directorio_raiz)
