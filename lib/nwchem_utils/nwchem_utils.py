import os

from installed_clients.KBaseReportClient import KBaseReport

# from lib.installed_clients.KBaseReportClient import KBaseReport

path = os.path.abspath(__file__)
template_folder = os.path.dirname(path)
template_folder = os.path.join(template_folder,'templates')
parameter_file = os.path.join(template_folder,'parameter_defaults.yaml')

def hello():
    print("hello")

    print("template_folder=", template_folder)

    print(os.listdir(template_folder))


def generate_report(params, callback_url):
    report_params = {'message': "Report will be built here",
                     'workspace_name': params.get('workspace_name')}

    kbase_report_client = KBaseReport(callback_url)
    report_info = kbase_report_client.create_extended_report(report_params)
    return report_info

#    def _mkdir_p(self, path):
#        """
#        _mkdir_p: make directory for given path
#        """
#        if not path:
#            return
#        try:
#            os.makedirs(path)
#        except OSError as exc:
#            if exc.errno == errno.EEXIST and os.path.isdir(path):
#                pass
#            else:
#                raise
