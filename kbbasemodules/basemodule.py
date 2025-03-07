from __future__ import absolute_import

import logging
import os
import copy
import json
import re
import time
import sys
import uuid
import requests
from os.path import exists
#from json import JSONEncoder
#class MyEncoder(JSONEncoder):
#def default_encoder(o):
#    return o.__dict__

logger = logging.getLogger(__name__)
logger.setLevel(
    logging.INFO
) 

class BaseModule:
    def __init__(self,name,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        #Initializing flexible container for client libraries which will be lazy loaded as needed
        self.version = "0.1.1.bm"
        self.clients = {}
        self.callback_url = callback
        for item in clients:
            self.clients[item] = clients[item]
        #Initializing config, name, and token
        self.config = config
        self.validate_args(self.config,[],{
            "max_retry":3,
            "workspace-url":"https://kbase.us/services/ws",
        })
        self.cached_to_obj_path = {}
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
    def initialize_call(self,method,params,print_params=False,no_print=[],no_prov_params=[]):
        if not self.initialized:
            self.obj_created = []
            self.input_objects = []
            self.method = method
            self.params = copy.deepcopy(params)
            for item in no_prov_params:
                if item in self.params:
                    del self.params[item]
            self.initialized = True
            if "workspace" in params:
                self.set_ws(params["workspace"])
            elif "output_workspace" in params:
                self.set_ws(params["output_workspace"])
            if print_params:
                saved_data = {}
                for item in no_print:
                    if item in params:
                        saved_data[item] = params[item]
                        del params[item]
                logger.info(method+":"+json.dumps(params,indent=4))
                for item in saved_data:
                    params[item] = saved_data[item]
                
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
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import Workspace
            else:
                from installed_clients.WorkspaceClient import Workspace
            self.clients["Workspace"] = Workspace(self.config["workspace-url"], token=self.token)
        return self.clients["Workspace"]
    
    def report_client(self):
        if "KBaseReport" not in self.clients:
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import KBaseReport
            else:
                from installed_clients.KBaseReportClient import KBaseReport
            self.clients["KBaseReport"] = KBaseReport(self.callback_url,token=self.token)
        return self.clients["KBaseReport"]
    
    def dfu_client(self):
        if "DataFileUtil" not in self.clients:
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import DataFileUtil
            else:
                from installed_clients.DataFileUtilClient import DataFileUtil
            self.clients["DataFileUtil"] = DataFileUtil(self.callback_url,token=self.token)
        return self.clients["DataFileUtil"]
    
    def gfu_client(self):
        if "GenomeFileUtil" not in self.clients:
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import GenomeFileUtil
            else:
                from installed_clients.GenomeFileUtilClient import GenomeFileUtil
            self.clients["GenomeFileUtil"] = GenomeFileUtil(self.callback_url,token=self.token)
        return self.clients["GenomeFileUtil"]
    
    def afu_client(self):
        if "AssemblyUtil" not in self.clients:
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import AssemblyUtil
            else:
                from installed_clients.AssemblyUtilClient import AssemblyUtil
            self.clients["AssemblyUtil"] = AssemblyUtil(self.callback_url,token=self.token)
        return self.clients["AssemblyUtil"]
    
    def rast_client(self):
        if "RAST_SDK" not in self.clients:
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import RAST_SDK
            else:
                from installed_clients.RAST_SDKClient import RAST_SDK
            self.clients["RAST_SDK"] = RAST_SDK(self.callback_url,token=self.token)
        return self.clients["RAST_SDK"]
    
    def anno_client(self,native_python_api=False):
        if "cb_annotation_ontology_api" not in self.clients:
            if native_python_api:
                from cb_annotation_ontology_api.annotation_ontology_api import AnnotationOntologyModule
                self.clients["cb_annotation_ontology_api"] = AnnotationOntologyModule("cb_annotation_ontology_api",{"data" :"/data/"},module_dir=self.module_dir+"/../cb_annotation_ontology_api",working_dir=self.working_dir,token=self.token,clients=self.clients,callback=self.callback_url)
            else:
                if "devenv" in self.config and self.config["devenv"] == "1":
                    from kbbasemodules import cb_annotation_ontology_api
                else:
                    from installed_clients.cb_annotation_ontology_apiClient import cb_annotation_ontology_api
                self.clients["cb_annotation_ontology_api"] = cb_annotation_ontology_api(self.callback_url,token=self.token)
        return self.clients["cb_annotation_ontology_api"]
    
    def handle_service(self):
        if "HandleService" not in self.clients:
            if "devenv" in self.config and self.config["devenv"] == "1":
                from kbbasemodules import AbstractHandle as HandleService
            else:
                from installed_clients.AbstractHandleClient import AbstractHandle as HandleService
            self.clients["HandleService"] = HandleService("https://kbase.us/services/handle_service", token=self.token)
        return self.clients["HandleService"]

    #########GENERAL UTILITY FUNCTIONS#######################
    def process_genome_list(self,input_references,workspace=None):
        ws_identities = []
        for ref in input_references:
            ws_identities.append(self.process_ws_ids(ref,workspace))
        output = self.ws_client().get_object_info3({"objects":ws_identities,"includeMetadata":1})
        output = output["infos"]
        output_references = []
        for info in output:
            if info[2].startswith("KBaseSearch.GenomeSet"):
                genomeset = self.get_object(self.wsinfo_to_ref(info))["data"]
                for label in genomeset["elements"]:
                    output_references.append(genomeset["elements"][label]["ref"])
            else:
                output_references.append(self.wsinfo_to_ref(info))
        return output_references
    
    def kb_version(self):
        wsclient = self.ws_client()
        if "appdev.kbase.us" in wsclient._client.url:
            return "dev"
        elif "/kbase.us" in wsclient._client.url:
            return "prod"
        elif "ci.kbase.us" in wsclient._client.url:
            return "ci"
        elif "next.kbase.us" in wsclient._client.url:
            return "next"
        else:
            return "unknown"
    
    def validate_args(self,params,required,defaults):
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

    def save_json(self,name,data):
        with open(self.working_dir+'/'+name+".json", 'w') as f:
            json.dump(data, f,indent=4,skipkeys=True)

    def load_json(self,name,default={}):
        if exists(self.working_dir+'/'+name+".json"):
            with open(self.working_dir+'/'+name+".json", 'r') as f:
                return json.load(f)
        return default
    
    def print_json_debug_file(self,filename,data):
        print("Printing debug file:",self.working_dir+'/'+filename)
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
            'service_ver': self.version
        }]
    
    def set_ws(self,workspace):
        if self.ws_id == workspace or self.ws_name == workspace:
            return 
        if not isinstance(workspace, str) or re.search('^\\d+$',workspace) != None:
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
        if len(id_or_ref.split(";")) > 1:
            objspec["to_obj_ref_path"] = id_or_ref.split(";")[0:-1]
            id_or_ref = id_or_ref.split(";")[-1]

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
            except Exception as e:
                logger.warning("Workspace get_objects2 call failed [%s:%s - %s]. Trying again! Error: %s", 
                       sys.exc_info()[0], sys.exc_info()[1], sys.exc_info()[2], str(e))
                tries += 1
                time.sleep(10)  # Give half second
        logger.warning("get_objects2 failed after multiple tries: %s", sys.exc_info()[0])
        raise

    def get_object_info(self, id_or_ref, ws=None):
        ws_identities = [self.process_ws_ids(id_or_ref, ws)]
        return self.ws_client().get_object_info(ws_identities,1)[0]

    def get_object(self, id_or_ref, ws=None):
        res = self.ws_get_objects({"objects": [self.process_ws_ids(id_or_ref, ws)]})
        if res is None:
            return None
        return res["data"][0]
    
    def save_genome_or_metagenome(self,objid,workspace,obj_json):
        self.set_ws(workspace)
        save_output = self.gfu_client().save_one_genome({
            "name" : objid,
            "data" : obj_json,
            "upgrade" : 1,
            "provenance" : self.provenance(),
            "hidden" : 0,
            "workspace" : self.ws_name
        });
        self.obj_created.append({"ref":self.create_ref(objid,self.ws_name),"description":""})
        return save_output["info"]
    
    def save_ws_object(self,objid,workspace,obj_json,obj_type):
        self.set_ws(workspace)
        params = {
            'id':self.ws_id,
            'objects': [{
                'data': obj_json,
                'name': objid,
                'type': obj_type,
                'meta': {},
                'provenance': self.provenance()
            }]
        }
        self.obj_created.append({"ref":self.create_ref(objid,self.ws_name),"description":""})
        return self.ws_client().save_objects(params)
    
    def wsinfo_to_ref(self,info):
        return str(info[6])+"/"+str(info[0])+"/"+str(info[4])
    
    def create_ref(self,id_or_ref,ws=None):
        if isinstance(id_or_ref, int):
            id_or_ref=str(id_or_ref)
        if len(id_or_ref.split("/")) > 1:
            return id_or_ref
        if isinstance(ws, int):
            ws=str(ws)
        return ws+"/"+id_or_ref
    
    def download_blob_file(self,handle_id,file_path,shock_url="https://kbase.us/services/shock-api"):
        headers = {'Authorization': 'OAuth ' + self.token}
        hs = self.handle_service()
        handles = hs.hids_to_handles([handle_id])
        shock_id = handles[0]['id']
        node_url = shock_url + '/node/' + shock_id
        r = requests.get(node_url, headers=headers, allow_redirects=True)
        errtxt = ('Error downloading file from shock ' +
                    'node {}: ').format(shock_id)
        if not r.ok:
            print(json.loads(r.content)['error'][0])
            return None
        resp_obj = r.json()
        size = resp_obj['data']['file']['size']
        if not size:
            print('Node {} has no file'.format(shock_id))
            return None
        node_file_name = resp_obj['data']['file']['name']
        attributes = resp_obj['data']['attributes']
        #Making the directory if it doesn't exist
        dir = os.path.dirname(file_path)
        os.makedirs(dir, exist_ok=True)
        #Adding filename to the end of the directory
        if os.path.isdir(file_path):
            file_path = os.path.join(file_path, node_file_name)
        with open(file_path, 'wb') as fhandle:
            with requests.get(node_url + '?download_raw', stream=True,
                                headers=headers, allow_redirects=True) as r:
                if not r.ok:
                    print(json.loads(r.content)['error'][0])
                    return None
                for chunk in r.iter_content(1024):
                    if not chunk:
                        break
                    fhandle.write(chunk)
        return file_path

    #########REPORT RELATED FUNCTIONS#######################
    def save_report_to_kbase(self,height=700,message="",warnings=[],file_links=[],summary_height=None):
        rootDir = self.working_dir+"/html/"
        files = [{'path': "/kb/module/work/tmp/html/",'name': "index.html",'description': 'HTML report'}]
        for dirName, subdirList, fileList in os.walk(rootDir):
            for fname in fileList:
                if fname != "index.html":
                    files.append({'path': dirName.replace(rootDir,"/kb/module/work/tmp/html/"),'name': fname,'description': 'Files related to HTML report'})
        report_name = self.method+"-"+str(uuid.uuid4())
        output = self.report_client().create_extended_report({
            'message': message,
            'warnings': warnings,            
            'html_links': files,
            'file_links': file_links,
            'direct_html_link_index': 0,
            'html_window_height': height,
            'objects_created': self.obj_created,
            'workspace_name': self.ws_name,
            'report_object_name': report_name,
            'summary_window_height': summary_height
        })
        return {"report_name":report_name,"report_ref":output["ref"],'workspace_name':self.ws_name}
    
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
                    "url": "data.json"
                },
                "columns": """+json.dumps(columns,indent=4)+"""
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