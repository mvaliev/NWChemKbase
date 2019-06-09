import io
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

print("path=",sys.path)
from decorators import memoize

def get_formula(mol):
    return Chem.rdMolDescriptors.CalcMolFormula(mol)

def get_mult(mol):
    return Descriptors.NumRadicalElectrons(mol) % 2 + 1

def get_charge(mol):
    return Chem.rdmolops.GetFormalCharge(mol)

def get_atom_symbols(mol):
    return [atom.GetSymbol() for atom in mol.GetAtoms()]

def get_atom_pos(mol,conf_id=0):

    c = mol.GetConformer(conf_id)
    return c.GetPositions()

def smiles_to_3Dmol(smiles):

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"failed creating molecule from smiles {smiles} \n")
        raise ValueError

    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol)

    AllChem.MMFFOptimizeMolecule(mol)

    return mol

def get_geom_block(mol, conf_id=0, filename=None):

    if filename is None:
        fp = io.StringIO()
    else:
        fp = open(filename, 'w', encoding="latin-1")

    name = get_atom_symbols(mol)
    pos =  get_atom_pos(mol,conf_id)

    for row, a in zip(pos, name):
        fp.write(a + '    ' + '  '.join('{: .12f}'.format(xyz) for xyz in row))
        fp.write('\n')

    try:
        output = fp.getvalue()
    except AttributeError:
        output = None

    fp.close()
    return output

def get_geom_xyz(mol,conf_id=0, filename=None):

    geom = get_geom_block(mol,conf_id=conf_id)

    natoms = mol.GetNumAtoms()
    formula = get_formula(mol)

    output = f'{natoms}\n{formula} conf-{conf_id}\n{geom}'

    if filename is not None:
        with open(filename, 'w', encoding="latin-1") as fp:
            fp.write(output)

    return output


if __name__ == '__main__':

    sm = '[OH2]'

    mol2 = smiles_to_3Dmol(sm)

    print(get_atom_pos(mol2))
    print(get_geom_block(mol2))
    print(get_geom_xyz(mol2))


