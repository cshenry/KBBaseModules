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
        self.version = "0.1.1.mm"
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.kbase_api = cobrakbase.KBaseAPI(token=token)
        #self.kbase_api = cobrakbase.KBaseCache(token=token,dev=True)
        self.kbase_api.ws_client = self.ws_client()
        #Loading default templates
        self.templates = {
            "core" : "NewKBaseModelTemplates/Core-V5.1",
            "gp" : "NewKBaseModelTemplates/GramPosModelTemplateV5",
            "gn" : "NewKBaseModelTemplates/GramNegModelTemplateV5",
            "ar" : "NewKBaseModelTemplates/ArchaeaTemplateV5",
            "grampos" : "NewKBaseModelTemplates/GramPosModelTemplateV5",
            "gramneg" : "NewKBaseModelTemplates/GramNegModelTemplateV5",
            "archaea" : "NewKBaseModelTemplates/ArchaeaTemplateV5",
            "old_grampos" : "NewKBaseModelTemplates/GramPosModelTemplateV3",
            "old_gramneg" : "NewKBaseModelTemplates/GramNegModelTemplateV3",
            "custom": None
        }
        
        #Setting ATP media
        if "ATP_media_workspace" not in self.config or not self.config["ATP_media_workspace"]:
            if self.kb_version() == "prod":
                self.config["ATP_media_workspace"] = "94026"
            elif self.kb_version() == "dev":
                self.config["ATP_media_workspace"] = "68393"
            else:
                logger.critical("KBase version not set up for modeling!")
    
    #################Utility functions#####################
    def process_media_list(self,media_list,default_media,workspace):
        #Retrieving media objects from references
        media_objects = []
        first = True
        #Cleaning out empty or invalid media references
        original_list = media_list
        media_list = []
        for media_ref in original_list:
            if len(media_ref) == 0:
                if first:
                    media_list.append(default_media)
                    first = False
                else:
                    print("Filtering out empty media reference")
            elif len(media_ref.split("/")) == 1:
                media_list.append(str(workspace)+"/"+media_ref)
            elif len(media_ref.split("/")) <= 3:
                media_list.append(media_ref)
            else:
                print(media_ref+" looks like an invalid workspace reference")
        #Making sure default gapfilling media is complete media
        if not media_list or len(media_list) == 0:
            media_list = [default_media]            
        #Retrieving media objects        
        for media_ref in media_list:  
            media = self.get_media(media_ref,None)
            media_objects.append(media)
        return media_objects
    
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
    
    def get_template(self,template_id,ws=None):
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
        mdlutl.wsid = objid
        #Saving attributes and getting model data
        if not isinstance(mdlutl.model,FBAModel):
            mdlutl.model = CobraModelConverter(mdlutl.model).build()
        mdlutl.save_attributes()
        data = mdlutl.model.get_data()
        #If the workspace is None, then saving data to file
        if not workspace:
            self.print_json_debug_file(mdlutl.wsid+".json",data)
        else:
            #Setting the workspace
            if workspace:
                self.set_ws(workspace)
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