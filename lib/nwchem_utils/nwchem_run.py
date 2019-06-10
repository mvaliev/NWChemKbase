import os
import subprocess


from nwchem_input import generate_input

from nwchem_mol import get_mult
from nwchem_mol import get_charge
from nwchem_mol import smiles_to_3Dmol
from nwchem_mol import get_geom_block

nwchem = ''
run_dir = '/Users/marat/tmp/nwchem_run'

try:
    nwchem = os.environ['NWCHEM_EXECUTABLE']
except KeyError:
    pass

try:
    run_dir = os.environ['NWCHEM_SIM_DIR']
except KeyError:
    pass

def generate_input_path():
    return os.path.join(run_dir,'nwchem-0.nw')

def generate_output_path():
    return os.path.join(run_dir,'nwchem-0.out')

def nwchem_run_smiles(smiles):

    output  = ''

    mol = smiles_to_3Dmol(smiles)

    geom = get_geom_block(mol)

    mult = get_mult(mol)
    charge = get_charge(mol)

    input_path = generate_input_path()
    output_path = generate_output_path()

    fo = open(output_path,'w',encoding="latin-1")
    print(input_path)
    generate_input(output,mult=mult,geometry=geom,charge=charge, filename=input_path,task='optimization')

    # completed = subprocess.run(['ls', run_dir])
    # print('returncode:', completed.returncode)
    # print(completed)

    completed = subprocess.run(['mkdir', os.path.join(run_dir,'perm')])
    print('completed:',completed)
    print('returncode:', completed.returncode)

    completed = subprocess.run(['mkdir', os.path.join(run_dir,'data')])
    print('completed:',completed)
    print('returncode:', completed.returncode)

    completed = subprocess.run(['mpirun', '-np', '2','--allow-run-as-root',nwchem,input_path],
                               encoding='latin-1',
                               # stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=run_dir,
                               stdout=fo, stderr=subprocess.STDOUT, cwd=run_dir,
                               shell=False, universal_newlines=True)
                               # stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=run_dir,
                               # shell=False, universal_newlines=True,text=True)
                               # shell=False)
    # completed = subprocess.run(['cat',input_path])
    fo.close()
    # print('completed:',completed)
    # print('returncode:', completed.returncode)
    # print('output:', completed.stdout)

    completed = subprocess.run(['cat',output_path],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
    # print('completed:',completed)
    # print('returncode:', completed.returncode)
    # print('output file:', completed.stdout)

    # completed = subprocess.run(['cat',input_path],
    #                            stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
    # print('completed:',completed)
    # print('returncode:', completed.returncode)
    # print('input file:', completed.stdout)


    # pattern = r'Total\s+\w+\s+[e,E]nergy\s+=\s+([-,+,\d,.]+)'
    # #
    # #
    # for line in completed.stdout.split('\n'):
    #     # print("line:",line)
    #     match = re.search(pattern, line)
    #     if match:
    #         print(match.groups())

    # pattern = r'Total\s+\w+\s+[e,E]nergy\s+=\s+([-,+,\d,.]+)'
    #
    # with open(output_path, 'r') as fp:
    #     for line in fp:
    #         # print(line)
    #         match = re.search(pattern, line)
    #         if match:
    #             print(match.groups())

    return input_path,output_path

if __name__ == '__main__':
    nwchem_run_smiles('[OH]')