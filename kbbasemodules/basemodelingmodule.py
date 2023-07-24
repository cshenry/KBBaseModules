from __future__ import absolute_import

import logging
import os
import sys
import json
import cobrakbase
from kbbasemodules.basemodule import BaseModule
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgrowthphenotypes import MSGrowthPhenotypes
from modelseedpy.core.msgenomeclassifier import MSGenomeClassifier
from cobrakbase.core.kbasefba.fbamodel_from_cobra import CobraModelConverter
from cobrakbase.core.kbasefba import FBAModel
from os.path import exists
import pickle

logger = logging.getLogger(__name__)

excluded_cpd = ["cpd22290","cpd11850"]

class BaseModelingModule(BaseModule):
    def __init__(self,name,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        BaseModule.__init__(self,name,config,module_dir,working_dir,token,clients,callback)
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.kbase_api = cobrakbase.KBaseAPI(token=token)
        #self.kbase_api = cobrakbase.KBaseCache(token=token,dev=True)
        self.kbase_api.ws_client = self.ws_client()
        #Loading default templates
        self.templates = {
            "core" : self.get_template("Core-V5.1","NewKBaseModelTemplates"),
            "gp" : None,
            "gn" : None,
            "custom": None
        }        
    
    #################Genome functions#####################
    def get_msgenome(self,id_or_ref,ws=None):
        genome = self.kbase_api.get_from_ws(id_or_ref,ws)
        self.input_objects.append(genome.info.reference)
        return genome
    
    def get_media(self,id_or_ref,ws=None):
        media = self.kbase_api.get_from_ws(id_or_ref,ws)
        media.id = media.info.id
        self.input_objects.append(media.info.reference)
        return media
    
    def get_phenotypeset(self,id_or_ref,ws=None,base_media=None, base_uptake=0, base_excretion=1000,global_atom_limits={}):
        kbphenoset = self.kbase_api.get_object(id_or_ref,ws)
        phenoset = MSGrowthPhenotypes.from_kbase_object(kbphenoset,self.kbase_api,base_media,base_uptake,base_excretion,global_atom_limits)
        return phenoset
    
    def get_model(self,id_or_ref,ws=None):
        mdlutl = MSModelUtil(self.kbase_api.get_from_ws(id_or_ref,ws))
        mdlutl.wsid = mdlutl.model.info.id
        #kbmodel = self.kbase_api.get_object(mdl_ref,None) #Should not have to do these three steps if the cobrakbase is working right
        #mdlutl.model.genome = self.kbase_api.get_from_ws(kbmodel["genome_ref"],None)
        #mdlutl.model.template = self.kbase_api.get_from_ws(kbmodel["template_ref"],None)
        self.input_objects.append(mdlutl.model.info.reference)
        return mdlutl
    
    #################Classifier functions#####################
    def get_classifier(self):
        cls_pickle = self.config["data"]+"/knn_ACNP_RAST_full_01_17_2023.pickle"
        cls_features = self.config["data"]+"/knn_ACNP_RAST_full_01_17_2023_features.json"
        #cls_pickle = self.module_dir+"/data/knn_ACNP_RAST_filter.pickle"
        #cls_features = self.module_dir+"/data/knn_ACNP_RAST_filter_features.json"
        with open(cls_pickle, 'rb') as fh:
            model_filter = pickle.load(fh)
        with open(cls_features, 'r') as fh:
            features = json.load(fh)
        return MSGenomeClassifier(model_filter, features)
    
    #################Template functions#####################
    def get_gs_template(self,template_id,ws,core_template):
        gs_template = self.get_template(template_id,ws)
        for cpd in core_template.compcompounds:
            if cpd.id not in gs_template.compcompounds:
                gs_template.compcompounds.append(cpd)
        for rxn in core_template.reactions:
            if rxn.id in gs_template.reactions:
                gs_template.reactions._replace_on_id(rxn)
            else:
                gs_template.reactions.append(rxn)
        for rxn in gs_template.reactions:
            for met in rxn.metabolites:
                if met.id[0:8] in excluded_cpd:
                    gs_template.reactions.remove(rxn)
        return gs_template
    
    def get_template(self,template_id,ws):
        template = self.kbase_api.get_from_ws(template_id,ws)
        #template = self.kbase_api.get_object(template_id,ws)
        #info = self.kbase_api.get_object_info(template_id,ws)
        #template = MSTemplateBuilder.from_dict(template).build()
        self.input_objects.append(template.info.reference)
        return template

    #################Save functions#####################
    def save_model(self,mdlutl,workspace=None,objid=None,suffix=None):
        #Setting the ID based on input
        if not suffix:
            suffix = ""
        if not objid:
            objid = mdlutl.wsid
        if not objid:
            logger.critical("Must provide an ID to save a model!")
        objid = objid+suffix
        
        #Setting the workspace
        if workspace:
            self.set_ws(workspace)
        
        #Saving final attributes and converting the model into KBase format if needed
        self.print_json_debug_file(objid+"-attributes.json",mdlutl.attributes)
        #if not isinstance(mdlutl.model,FBAModel):
            #mdlutl.model = CobraModelConverter(mdlutl.model,mdlutl.model.genome, mdlutl.model.template).build()
        mdlutl.save_attributes()
        data = mdlutl.model.get_data()
        
        #Setting provenance and saving model using workspace API
        mdlutl.create_kb_gapfilling_data(data,self.config["ATP_media_workspace"])
        params = {
            'id':self.ws_id,
            'objects': [{
                'data': data,
                'name': objid,
                'type': "KBaseFBA.FBAModel",
                'meta': {},
                'provenance': self.provenance()
            }]
        }
        self.ws_client().save_objects(params)
        self.obj_created.append({"ref":self.create_ref(objid,self.ws_name),"description":""})