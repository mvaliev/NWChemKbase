#!/usr/bin/env python

from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import argparse


parser = argparse.ArgumentParser()

parser.add_argument("input")

args = parser.parse_args()

input = args.input

mol = Chem.MolFromSmiles(input)

molh=Chem.AddHs(mol)

AllChem.EmbedMolecule(molh)

AllChem.MMFFOptimizeMolecule(molh)

pdb = Chem.MolToPDBBlock(molh) 
charge=Chem.rdmolops.GetFormalCharge(molh)
spin=int(Descriptors.NumRadicalElectrons(molh))+1

with open('par.dat','w') as f:
  f.write('charge=%d\n'%charge)
  f.write('mult=%d\n'%spin)


with open('mol.pdb','w') as f:
  f.write(pdb)
