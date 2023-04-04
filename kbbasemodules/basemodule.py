from __future__ import absolute_import

import logging
import os
import copy
import json
import re
import time
import sys
import uuid
from os.path import exists
#from json import JSONEncoder
#class MyEncoder(JSONEncoder):
#def default_encoder(o):
#    return o.__dict__

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
    level=logging.INFO)

class BaseModule:
    def __init__(self,name,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        #Initializing flexible container for client libraries which will be lazy loaded as needed
        self.clients = {}
        self.callback_url = callback
        for item in clients:
            self.clients[item] = clients[item]
        #Initializing config, name, and token
        self.config = config
        self.validate_args(self.config,[],{
            "max_retry":3,
            "version":0,
            "workspace-url":"https://kbase.us/services/ws",
        })
        self.token = token
        self.name = name
        self.module_dir = module_dir
        #Initializing working directory if specified, otherwise using config scratch
        if working_dir:
            self.working_dir = working_dir
        else:
            self.working_dir = config['scratch']
        self.reset_attributes()
    
    #########METHOD CALL INITIALIZATION FUNCTIONS#######################
    def initialize_call(self,method,params,print_params=False):
        if not self.initialized:
            self.obj_created = []
            self.input_objects = []
            self.method = method
            self.params = copy.deepcopy(params)
            self.initialized = True
            if "workspace" in params:
                self.set_ws(params["workspace"])
            if print_params:
                logger.info(method+":"+json.dumps(params,indent=4))
                
    def reset_attributes(self):
        #Initializing stores tracking objects created and input objects
        self.obj_created = []
        self.input_objects = []
        #Initializing attributes tracking method data to support provencance and context
        self.method = None
        self.params = {}
        self.initialized = False
        self.ws_id = None
        self.ws_name = None
        #Computing timestamp
        ts = time.gmtime()
        self.timestamp = time.strftime("%Y-%m-%d %H:%M:%S", ts)
    
    #########CLIENT RETRIEVAL AND INITIALIZATION FUNCTIONS#######################
    def ws_client(self):
        if "Workspace" not in self.clients:
            from installed_clients.WorkspaceClient import Workspace
            self.clients["Workspace"] = Workspace(self.config["workspace-url"], token=self.token)
        return self.clients["Workspace"]
    
    def report_client(self):
        if "KBaseReport" not in self.clients:
            from installed_clients.KBaseReportClient import KBaseReport
            self.clients["KBaseReport"] = KBaseReport(self.callback_url,token=self.token)
        return self.clients["KBaseReport"]
    
    def dfu_client(self):
        if "DataFileUtil" not in self.clients:
            from installed_clients.DataFileUtilClient import DataFileUtil
            self.clients["DataFileUtil"] = DataFileUtil(self.callback_url,token=self.token)
        return self.clients["DataFileUtil"]

    #########GENERAL UTILITY FUNCTIONS#######################
    def validate_args(self,params,required,defaults):
        #print("One:"+json.dumps(params,indent=4)+"\n\n")
        for item in required:
            if item not in params:
                raise ValueError('Required argument '+item+' is missing!')
        for key in defaults:
            if key not in params:
                params[key] = defaults[key]
        return params
    
    def transfer_outputs(self,output,api_output,key_list):
        for key in key_list:
            if key in api_output:
                output[key] = api_output[key]
                
    def print_json_debug_file(self,filename,data):
        with open(self.working_dir+'/'+filename, 'w') as f:
            json.dump(data, f,indent=4,skipkeys=True)
    
    #########WORKSPACE RELATED FUNCTIONS#######################
    def provenance(self):
        return [{
            'description': self.method,
            'input_ws_objects': self.input_objects,
            'method': self.method,
            'script_command_line': "",
            'method_params': [self.params],
            'service': self.name,
            'service_ver': self.config["version"],
            # 'time': '2015-12-15T22:58:55+0000'
        }]
    
    def set_ws(self,workspace):
        if self.ws_id == workspace or self.ws_name == workspace:
            return 
        if not isinstance(workspace, str) or re.search('^\d+$',workspace) != None:
            if isinstance(workspace, str):
                workspace = int(workspace)
            self.ws_id = workspace
            info = self.ws_client().get_workspace_info({"id":workspace})
            self.ws_name = info[1]
        else:
            self.ws_name = workspace
            info = self.ws_client().get_workspace_info({"workspace":workspace})
            self.ws_id = info[0]
    
    def process_ws_ids(self,id_or_ref,workspace=None,no_ref=False):
        """
        IDs should always be processed through this function so we can interchangeably use
        refs, IDs, and names for workspaces and objects
        """
        objspec = {}
        if len(id_or_ref.split("/")) > 1:
            if no_ref:
                array = id_or_ref.split("/")
                workspace = array[0]
                id_or_ref = array[1]
            else:
                objspec["ref"] = id_or_ref
                
        if "ref" not in objspec:
            if isinstance(workspace, int):
                objspec['wsid'] = workspace
            else:
                objspec['workspace'] = workspace
            if isinstance(id_or_ref, int):
                objspec['objid'] = id_or_ref
            else:
                objspec['name'] = id_or_ref
        return objspec
          
    def ws_get_objects(self, args):
        """
        All functions calling get_objects2 should call this function to ensure they get the retry
        code because workspace periodically times out
        :param args:
        :return:
        """
        tries = 0
        while tries < self.config["max_retry"]:
            try:
                return self.ws_client().get_objects2(args)
            except:
                logger.warning("Workspace get_objects2 call failed [%s:%s - %s]. Trying again!")
                tries += 1
                time.sleep(500)  # Give half second
        logger.warning("get_objects2 failed after multiple tries: %s", sys.exc_info()[0])
        raise

    def get_object(self, id_or_ref, ws=None):
        res = self.ws_get_objects({"objects": [self.process_ws_ids(id_or_ref, ws)]})
        if res is None:
            return None
        return res["data"][0]
    
    def create_ref(self,id_or_ref,ws=None):
        if isinstance(id_or_ref, int):
            id_or_ref=str(id_or_ref)
        if len(id_or_ref.split("/")) > 1:
            return id_or_ref
        if isinstance(ws, int):
            ws=str(ws)
        return ws+"/"+id_or_ref
    
    #########REPORT RELATED FUNCTIONS#######################
    def save_report_to_kbase(self,height=700,message=""):
        files = []
        rootDir = self.working_dir+"/html/"
        for dirName, subdirList, fileList in os.walk(rootDir):
            print('Found directory: %s' % dirName)
            for fname in fileList:
                print('\t%s' % fname)
                files.append({'path': dirName,'name': fname,'description': 'HTML report for PDB upload'})
        report_name = self.method+"-"+str(uuid.uuid4())
        output = self.report_client().create_extended_report({'message': message,
                         'html_links': files,
                         'direct_html_link_index': 0,
                         'html_window_height': height,
                         'objects_created': self.obj_created,
                         'workspace_name': self.ws_name,
                         'report_object_name': report_name})
        return {"report_name":report_name,"report_ref":output["ref"],'workspace_name':self.api.ws_name}
    
    def build_dataframe_report(self,table,column_list):        
        #Convert columns to this format:
        columns = []
        for item in column_list:
            columns.append({"data":item})
        #for index, row in table.iterrows():
        #    pass
        json_str = table.to_json(orient='records')
        #columns=column_list
        html_data = """
    <html>
    <header>
        <link href="https://cdn.datatables.net/1.11.5/css/jquery.dataTables.min.css" rel="stylesheet">
    </header>
    <body>
    <script src="https://code.jquery.com/jquery-3.6.0.slim.min.js" integrity="sha256-u7e5khyithlIdTpu22PHhENmPcRdFiHRjhAuHcs05RI=" crossorigin="anonymous"></script>
    <script type="text/javascript" src="https://cdn.datatables.net/1.11.5/js/jquery.dataTables.min.js"></script>
    <script>
        $(document).ready(function() {
            $('#example').DataTable( {
                "ajax": {
                    "url": "data.json",
                    "dataSrc": ""
                },
                "columns": {json.dumps(columns,indent=4)}
            } );
        } );
    </script>
    </body>
    </html>
    """
        os.makedirs(self.working_dir+"/html", exist_ok=True)
        with open(self.working_dir+"/html/index.html", 'w') as f:
            f.write(html_data)
        with open(self.working_dir+"/html/data.json", 'w') as f:
            f.write(json_str)