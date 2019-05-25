import sys

import os

import yaml
from jinja2 import Environment, FileSystemLoader, StrictUndefined, UndefinedError, TemplateNotFound

print(sys.path)

path = os.path.abspath(__file__)
template_folder = os.path.dirname(path)
template_folder = os.path.join(template_folder, 'templates')
parameter_file = os.path.join(template_folder, 'parameter_defaults.yaml')


def generate_input(output, filename=None, **kwargs):
    '''
    generate nwchem input file
    '''

    with open(parameter_file, 'r') as stream:
        params = yaml.safe_load(stream)

    params.update(kwargs)

    templates_dir = [template_folder]

    env = Environment(loader=FileSystemLoader(templates_dir), trim_blocks=False,
                      # undefined=StrictUndefined,
                      lstrip_blocks=True, line_statement_prefix="%",
                      line_comment_prefix="%%", extensions=['jinja2.ext.do'])

    template_file = params['task'] + '.nwt'

    try:
        template = env.get_template(template_file)
    except TemplateNotFound as e:
        print('cannot match task %s to available templates' % params['task'])
        raise

    buffer = template.render(params)

    if filename is not None:
        with open(filename, 'w', encoding="latin-1")  as fp:
            fp.write(buffer)
    else:
        return buffer


if __name__ == '__main__':
    from nwchem_mol import get_mult
    from nwchem_mol import get_charge
    from nwchem_mol import smiles_to_3Dmol
    from nwchem_mol import get_geom_block

    output = ''
    # generate_input(output)
    sm = '[OH]'

    mol = smiles_to_3Dmol(sm)

    geom = get_geom_block(mol)

    mult = get_mult(mol)
    charge = get_charge(mol)

    generate_input(output, mult=mult, geometry=geom, charge=charge, filename='energy.nw')
