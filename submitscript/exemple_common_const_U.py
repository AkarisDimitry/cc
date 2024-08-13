from common import const_U_relax
from ase.calculators.vasp import Vasp
from ase.io import read
import os,sys,pickle,math

kpts = (2,2,1)
pot_des = 1.230
ediffg = 0.05

try:
    atoms = read('CONTCAR')
except:
    atoms = read('init.traj')

calc = Vasp(txt='-',
            encut=500,
            xc='PBE',
            gga='RP',
            kpts  = kpts,
            ncore=16,
            kpar=1,
            gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
            ismear=0,
            algo = 'Normal',
            nelm=250,
            sigma = 0.20,
            lorbit=11,
            ibrion=2,
            ediffg=-1*ediffg,  # forces
            ediff=1e-4,
            prec='Accurate',
            nsw=1000,
            lcharg=False,
            lsol=True,
            lambda_d_k=3.0,
            tau=0.0,
            lasph=True)

const_U_relax(atoms=atoms,calc=calc,desired_U=pot_des)
