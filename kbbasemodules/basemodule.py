from __future__ import absolute_import

import logging
import os
import copy
import json
import re
import time
import sys
from os.path import exists

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
    level=logging.INFO)

class BaseModule:
    def __init__(self,name,ws_client,working_dir,config):
        self.ws_client = ws_client
        self.working_dir = working_dir
        self.obj_created = []
        self.input_objects = []
        self.method = None
        self.params = None
        self.name = name
        self.config = config
        self.initialized = False
        self.ws_id = None
        self.ws_name = None
        self.timestamp = time.time()
        self.validate_args(self.config,[],{
            "max_retry":3,
            "version":0
        })
    
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
                print(json.dumps(params,indent=4))
    
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
            json.dump(data, f)
    
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
            info = self.ws_client.get_workspace_info({"id":workspace})
            self.ws_name = info[1]
        else:
            self.ws_name = workspace
            info = self.ws_client.get_workspace_info({"workspace":workspace})
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
                return self.ws_client.get_objects2(args)
            except:
                logger.warning("Workspace get_objects2 call failed [%s:%s - %s]. Trying again!")
                tries += 1
                time.sleep(500)  # Give half second
        logger.warning("get_objects2 failed after multiple tries: %s", sys.exc_info()[0])
        raise

    def get_object(self, id_or_ref, ws=None):
        res = self.ws_get_objects({"objects": [self.process_ws_ids(id_or_ref, ws)]})
        self.print_json_debug_file(res,"TestGenome.json")
        if res is None:
            return None
        return res["data"][0]