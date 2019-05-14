# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import subprocess as _subprocess
import uuid
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

from installed_clients.KBaseReportClient import KBaseReport
from nwchem_utils.nwchem_utils import hello, generate_report
#END_HEADER


class NWChemKbase:
    '''
    Module Name:
    NWChemKbase

    Module Description:
    A KBase module: NWChemKbase
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/mvaliev/NWChemKbase.git"
    GIT_COMMIT_HASH = "eb1d1a8b6c05965372b248c5c0f0940c53cc687b"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.shared_folder = config['scratch']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.scratch = config['scratch']
        print('IN INIT DIRECTORY')
        #END_CONSTRUCTOR
        pass


    def run_NWChemKbase(self, ctx, params):
        """
        This example function accepts any number of parameters and returns results in a KBaseReport
        :param params: instance of mapping from String to unspecified object
        :returns: instance of type "ReportResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_NWChemKbase

        result_directory = os.path.join(self.scratch, str(uuid.uuid4()))
#        self._mkdir_p(result_directory)

        #s = _subprocess.call(["echo $NWCHEM_EXECUTABLE"])
        s = _subprocess.call(["ls"])
        print(s)
        nwchem = os.environ['NWCHEM_EXECUTABLE']
        run_nwchem = os.environ['NWCHEM_BIN']+"/run_nwchem"

        print (run_nwchem)
        nwchem_output = os.environ['NWCHEM_SIM_DIR']+"/nwchem.out"
        s = _subprocess.call([run_nwchem,params['smiles_string']])
        print("output file ",nwchem_output)
        output = _subprocess.run(['cat',nwchem_output],stdout=_subprocess.PIPE)
        energy = _subprocess.run(['grep','Total DFT',nwchem_output],stdout=_subprocess.PIPE)
        energy = energy.stdout.decode('utf-8')
        #print("nwchem output=",output.stdout.decode('utf-8'))
        print("params=",params)
        mol = Chem.MolFromSmiles(params['smiles_string'])
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
        print("formula",formula)
        text_message = "".join([
           'NWChem single point energy calculation:\n',
           'Molecular SMILE string:',
           str(params['smiles_string']),
           '\n',
           'Molecular Formula: ',formula,
           '\n', str(energy).strip(),'\n'
        ])
        print("\n------------------")
        print(text_message)
        hello()
        print("------------------\n")


#        reportObj = {
#            'objects_created': [],
#            'text_message': text_message
#        }
#        report = KBaseReport(self.callback_url)
#        report_info = report.create({'report': reportObj, 'workspace_name': params['workspace_name']})
        # report_params = {'message': text_message,
        #                  'workspace_name': params.get('workspace_name')}

        # kbase_report_client = KBaseReport(self.callback_url)
        # report_info = kbase_report_client.create_extended_report(report_params)

        report_info = generate_report(params,self.callback_url )

        # STEP 6: contruct the output to send back
        output = {'report_name': report_info['name'],
                  'report_ref': report_info['ref'],
                  }

#        report = KBaseReport(self.callback_url)
#        report_info = report.create_extended_report({'report': {'objects_created':[],
#                                                'message': text_message},
#                                                'workspace_name': params['workspace_name']})
#        output = {
#            'report_name': report_info['name'],
#            'report_ref': report_info['ref'],
#        }

        #END run_NWChemKbase

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_NWChemKbase return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
