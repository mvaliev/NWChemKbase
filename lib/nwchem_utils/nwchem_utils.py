
from installed_clients.KBaseReportClient import KBaseReport

def hello():

   print("hello")


def generate_report(params,callback_url):

   report_params = {'message': "Just testing here",
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
