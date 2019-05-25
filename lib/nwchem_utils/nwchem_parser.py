import collections
import re
import sys
import os
import io


import rdkit

from nwchem_mol import get_mult

path = os.path.abspath(__file__)
data_folder = os.path.dirname(path)
data_folder = os.path.join(data_folder,'data')

Atom = collections.namedtuple('Atom', 'id tag x y z')

pattern = {
    'header': r'( +Northwest\W+[\w,\W]+99352)',
    'general': r'General\s+Information\W+-([\w,\W]+?Spin[^\n]+)',
    'citation': r'(^ +CITA[\w,\W]+)AU'
}

def get_chunk(output_path,pattern,center=False,shift=False, last_line=-1):

    with open(output_path, 'r') as fp:
        buffer = fp.read()

    # match = re.findall(pattern, buffer, re.MULTILINE | re.IGNORECASE)
    match = re.findall(pattern, buffer, re.IGNORECASE)

    # center = False
    # shift = True
    buffer = match[-1].split('\n')
    output = ''
    for line in buffer[:last_line]:
        if center:
            line = line.strip()
            factor = 30 - int(len(line)/2)
            line = ' '*factor+ line.strip()
        elif shift:
            line = line.strip()
            factor = 0
            line = ' '*factor+ line.strip()
        output = output + f'{line}\n'

    return output


def general_information(output_path):
    # pattern = r'General\s+Information\W+([\w,\W]+Spin[^\n]+)'
    #
    # with open(output_path, 'r') as fp:
    #     str = fp.read()
    #     # print(str)
    #     match = re.search(pattern,str,re.MULTILINE)
    #     print(match.end())

    pattern = r'General\s+Information\W+-([\w,\W]+?Spin[^\n]+)'

    pattern = r'( +Northwest\W+[\w,\W]+99352)'
    with open(output_path, 'r') as fp:
        buffer = fp.read()

    # match = re.findall(pattern, buffer, re.MULTILINE | re.IGNORECASE)
    match = re.findall(pattern, buffer, re.IGNORECASE)


    print(match[-1])

    print(len(match))


def parse_total_energy(output_path):

    pattern = r'Total\s+\w+\s+[e,E]nergy\s+=\s+([-,+,\d,.]+)'

    energy = []
    with open(output_path, 'r') as fp:
        for line in fp:
            # print(line)
            m = re.search(pattern, line)
            if m:
                energy.append(m.group(1))

    print(energy)

def parse_task(output_path):
    pattern = r'task\s+(?P<theory>\w+)\s+(?P<task>\w+)'
    with open(output_path, 'r') as fp:
        for line in fp:
            # print(line)
            match = re.search(pattern, line)
            if match:
                return match.groups()


def parse_summary():
    pass

def extract_geometry(output_path,image=None):
    fp = open(output_path, 'r')

    buffer = fp.read()

    pattern = r'Output\s+coordinates.*\s+\n\s+No.\s+Tag.*\n\s*-.*((?:\n\s+\d+.*)+)\s+Atomic'
    match = re.findall(pattern, buffer)
    geom_list = []
    for g in match:
        geom_block = g.strip().split('\n')
        geom = []
        for id, xyz_line in enumerate(geom_block,start=1):
            _,tag,_,x,y,z = xyz_line.split()
            geom.append(Atom(id, tag, float(x), float(y), float(z)))
            # print(tokens)
            # tag = tokens[1]
            # print(tag)
            # coords   = [float(a) for a in tokens[3:6]]
            # print(coords)
            # tag_coords = [tag] + coords
            # print(tag_coords)
        geom_list.append(geom)

    # for a in geom_list[-1]:
    #     print(a)

    fmt_string = "{record:<6s}{a.id:>5d}{blank:>1s}{a.tag:<4s}{altloc:>1s}" \
                 "{groupname:>3s}{blank:>1s}{chainid:>1s}{groupid:>4}" \
                 "{blank:>4s}" \
                 "{a.x:>8.3f}{a.y:>8.3f}{a.z:>8.3f}" \
                 "{occupancy:>6s}{tempfactor:>6s}" \
                 "{blank:>10s}{a.tag:>2s}\n"

    fp = io.StringIO()

    for a in geom_list[-1]:
        fp.write(fmt_string.format(a=a,
                             id=1,
                             groupname='UNK',
                             groupid=1,
                             record="ATOM",
                             blank="",
                             altloc="",
                             chainid="",
                             occupancy="",
                             tempfactor=""))


    return fp.getvalue()
    # return fmt_string.format(a=a,
    #                          id=id or a.id,
    #                          groupname=groupname or a.groupname,
    #                          groupid=groupid or a.groupid,
    #                          record="ATOM",
    #                          blank="",
    #                          altloc="",
    #                          chainid="",
    #                          occupancy="",
    #                          tempfactor="")


def generate_summary(output_file):

    header = get_chunk(output_file,pattern['header'],last_line=3,shift=True)

    theory, task = parse_task(output_file)
    if 'opt' in task:
        calc_type = 'OPTIMIZATION CALCULATION'
    if 'energy' in task:
        calc_type = 'ENERGY CALCULATION'

    general = get_chunk(output_file,pattern['general'],shift=True)

    return f'{header}\n{calc_type}\n{general}'


if __name__ == '__main__':

    from rdkit.Chem import Descriptors
    from rdkit.Chem.rdmolfiles import MolFromPDBBlock
    from rdkit.Chem import Draw
    output = os.path.join(data_folder,'dft-optimization.out')
    print(generate_summary(output))
    # # general_information(output)
    # print(get_chunk(output,pattern['header'],last_line=3,shift=True))
    #
    # theory, task = parse_task(output)
    # if 'opt' in task:
    #     print(f'OPTIMIZATION CALCULATION')
    #
    # print(get_chunk(output,pattern['general'],shift=True))
    #
    #
    #
    # # parse_total_energy(output)
    #
    #
    # # buffer = extract_geometry(output)
    # #
    # # print(buffer)
    # # mol = MolFromPDBBlock(buffer)
    # # print(mol)
    # #
    # #
    # # # def get_mult(mol):
    # # #     try:
    # # #         mult = mol.mult
    # # #     except AttributeError:
    # # #         print("computing")
    # # #         mult = Descriptors.NumRadicalElectrons(mol) % 2 + 1
    # # #         mol.mult = mult
    # # #
    # # #     return mult
    # #
    # # print('mult=',get_mult(mol))
    # # print('mult=',get_mult(mol))
    # # print('mult=', get_mult(mol))
    # # # Draw.MolToFile(mol, 'mol.png')
