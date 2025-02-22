from __future__ import absolute_import

import logging
import os
import sys
import json
import cobrakbase
from kbbasemodules.basemodule import BaseModule
import cobra
from cobra import Model, Reaction, Metabolite
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.core.msgenome import MSGenome
from modelseedpy.core.msfba import MSFBA
from modelseedpy.biochem.modelseed_biochem import ModelSEEDBiochem
from modelseedpy.biochem.modelseed_reaction import ModelSEEDReaction2
from modelseedpy.biochem.modelseed_compound import ModelSEEDCompound2
from modelseedpy.core.annotationontology import AnnotationOntology
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
        BaseModule.__init__(self,name,config,module_dir=module_dir,working_dir=working_dir,token=token,clients=clients,callback=callback)
        self.version = "0.1.1.mm"
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        self.kbase_api = cobrakbase.KBaseAPI(token=token)
        #self.kbase_api = cobrakbase.KBaseCache(token=token,dev=True)
        self.kbase_api.ws_client = self.ws_client()
        #Loading default templates
        self.templates = {
            "core" : "NewKBaseModelTemplates/Core-V5.2",
            "gp" : "NewKBaseModelTemplates/GramPosModelTemplateV6",
            "gn" : "NewKBaseModelTemplates/GramNegModelTemplateV6",
            "ar" : "NewKBaseModelTemplates/ArchaeaTemplateV6",
            "grampos" : "NewKBaseModelTemplates/GramPosModelTemplateV6",
            "gramneg" : "NewKBaseModelTemplates/GramNegModelTemplateV6",
            "archaea" : "NewKBaseModelTemplates/ArchaeaTemplateV6",
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
        if not media_list:
            media_list = []
        #Retrieving media objects from references
        media_objects = []
        first = True
        #Cleaning out empty or invalid media references
        original_list = media_list
        media_list = []
        for media_ref in original_list:
            if not media_ref or len(media_ref) == 0:
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
    def annotate_genome_with_rast(self,genome_id,ws=None,output_ws=None):
        if not output_ws:
            output_ws = ws
        rast_client = self.rast_client()
        output = rast_client.annotate_genome({
            "workspace": output_ws,
            "input_genome":genome_id,
            "output_genome":genome_id+".RAST"
        })
        return output["workspace"]+"/"+output["id"]
    
    def get_msgenome_from_ontology(self,id_or_ref,ws=None,native_python_api=False,output_ws=None):
        annoapi = self.anno_client(native_python_api=native_python_api)
        gen_ref = self.create_ref(id_or_ref,ws)
        genome_info = self.get_object_info(gen_ref)
        annoont = AnnotationOntology.from_kbase_data(annoapi.get_annotation_ontology_events({
            "input_ref" : gen_ref
        }),gen_ref,self.module_dir+"/data/")
        gene_term_hash = annoont.get_gene_term_hash(ontologies=["SSO"])
        if len(gene_term_hash) == 0:
            logger.warning("Genome has not been annotated with RAST! Reannotating genome with RAST!")
            gen_ref = self.annotate_genome_with_rast(genome_info[1],genome_info[6],output_ws)
            annoont = AnnotationOntology.from_kbase_data(annoapi.get_annotation_ontology_events({
                "input_ref" : gen_ref
            }),gen_ref,self.module_dir+"/data/")
        annoont.info = genome_info
        wsgenome = self.get_msgenome(gen_ref,ws)
        genome = annoont.get_msgenome()
        for ftr in wsgenome.features:
            for func in ftr.functions:
                if ftr.id in genome.features:
                    genome.features.get_by_id(ftr.id).add_ontology_term("RAST",func)
                else:
                    newftr = genome.create_new_feature(ftr.id,"")
                    newftr.add_ontology_term("RAST",func)
        genome.id = genome_info[1]
        genome.scientific_name = genome_info[10]["Name"]
        return genome

    def get_expression_objs(self,expression_refs,genome_objs):
        genomes_to_models_hash = {}
        for mdl in genome_objs:
            genomes_to_models_hash[genome_objs[mdl]] = mdl
        ftrhash = {}
        expression_objs = {}
        for genome_obj in genomes_to_models_hash:
            for ftr in genome_obj.features:
                ftrhash[ftr.id] = genome_obj
        for expression_ref in expression_refs:
            expression_obj = self.kbase_api.get_from_ws(expression_ref,None)
            row_ids = expression_obj.row_ids
            genome_obj_count = {}
            for ftr_id in row_ids:
                if ftr_id in ftrhash:
                    if ftrhash[ftr_id] not in genome_obj_count:
                        genome_obj_count[ftrhash[ftr_id]] = 0
                    genome_obj_count[ftrhash[ftr_id]] += 1
            best_count = None
            best_genome = None
            for genome_obj in genome_obj_count:
                if best_genome == None or genome_obj_count[genome_obj] > best_count:
                    best_genome = genome_obj
                    best_count = genome_obj_count[genome_obj]
            if best_genome:
                expression_objs[genomes_to_models_hash[best_genome]] = expression_obj.data    
        return expression_objs
    
    def column_ids(self):
        return self.data.columns.to_list()

    @property
    def row_ids(self):
        return self.data.index.to_list()

    def get_msgenome(self,id_or_ref,ws=None):
        genome = self.kbase_api.get_from_ws(id_or_ref,ws)
        genome.id = genome.info.id
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
    
    def get_model(self,id_or_ref,ws=None,is_json_file=False):
        if is_json_file:
            return MSModelUtil.build_from_kbase_json_file(id_or_ref)
        mdlutl = MSModelUtil(self.kbase_api.get_from_ws(id_or_ref,ws))
        mdlutl.wsid = mdlutl.model.info.id
        self.input_objects.append(mdlutl.model.info.reference)
        return mdlutl
    
    def extend_model_with_other_ontologies(self,mdlutl,anno_ont,builder,prioritized_event_list=None,ontologies=None,merge_all=True):
        gene_term_hash = anno_ont.get_gene_term_hash(
            prioritized_event_list, ontologies, merge_all, False
        )
        residual_reaction_gene_hash = {}
        for gene in gene_term_hash:
            for term in gene_term_hash[gene]:
                if term.ontology.id != "SSO":
                    for rxn_id in term.msrxns:
                        if rxn_id not in residual_reaction_gene_hash:
                            residual_reaction_gene_hash[rxn_id] = {}
                        if gene not in residual_reaction_gene_hash[rxn_id]:
                            residual_reaction_gene_hash[rxn_id][gene] = []
                        residual_reaction_gene_hash[rxn_id][gene] = gene_term_hash[
                            gene
                        ][term]

        reactions = []
        SBO_ANNOTATION = "sbo"
        modelseeddb = ModelSEEDBiochem.get()
        biochemdbrxn = False
        for rxn_id in residual_reaction_gene_hash:
            if rxn_id + "_c0" not in mdlutl.model.reactions:
                reaction = None
                template_reaction = None
                if rxn_id + "_c" in mdlutl.model.template.reactions:
                    template_reaction = mdlutl.model.template.reactions.get_by_id(rxn_id + "_c")
                elif rxn_id in modelseeddb.reactions:
                    rxnobj = modelseeddb.reactions.get_by_id(rxn_id)
                    if "MI" not in rxnobj.status and "CI" not in rxnobj.status:
                        #mdlutl.add_ms_reaction({rxn_id:"c0"}, compartment_trans=["c0", "e0"])
                        template_reaction = rxnobj.to_template_reaction({0: "c", 1: "e"})
                        biochemdbrxn = True
                if template_reaction:
                    for m in template_reaction.metabolites:
                        if m.compartment not in builder.compartments:
                            builder.compartments[
                                m.compartment
                            ] = builder.template.compartments.get_by_id(m.compartment)
                        if m.id not in builder.template_species_to_model_species:
                            model_metabolite = m.to_metabolite(builder.index)
                            builder.template_species_to_model_species[
                                m.id
                            ] = model_metabolite
                            builder.base_model.add_metabolites([model_metabolite])
                    if biochemdbrxn:
                        pass
                        #template_reaction.add_metabolites({})
                    reaction = template_reaction.to_reaction(
                        builder.base_model, builder.index
                    )
                    gpr = ""
                    probability = None
                    for gene in residual_reaction_gene_hash[rxn_id]:
                        for item in residual_reaction_gene_hash[rxn_id][gene]:
                            if "scores" in item:
                                if "probability" in item["scores"]:
                                    if (
                                        not probability
                                        or item["scores"]["probability"] > probability
                                    ):
                                        probability = item["scores"]["probability"]
                        if len(gpr) > 0:
                            gpr += " or "
                        gpr += gene.id
                    if probability != None and hasattr(reaction, "probability"):
                        reaction.probability = probability
                    reaction.gene_reaction_rule = gpr
                    reaction.annotation[SBO_ANNOTATION] = "SBO:0000176"
                    reactions.append(reaction)
                if not reaction:
                    print("Reaction ", rxn_id, " not found in template or database!")
            else:
                rxn = mdlutl.model.reactions.get_by_id(rxn_id + "_c0")
                gpr = rxn.gene_reaction_rule
                probability = None
                for gene in residual_reaction_gene_hash[rxn_id]:
                    for item in residual_reaction_gene_hash[rxn_id][gene]:
                        if "scores" in item:
                            if "probability" in item["scores"]:
                                if (
                                    not probability
                                    or item["scores"]["probability"] > probability
                                ):
                                    probability = item["scores"]["probability"]
                    if len(gpr) > 0:
                        gpr += " or "
                    gpr += gene.id
                if probability != None and hasattr(rxn, "probability"):
                    rxn.probability = probability
                rxn.gene_reaction_rule = gpr
        mdlutl.model.add_reactions(reactions)
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
    
    def save_phenotypeset(self,data,workspace,objid):
        self.set_ws(workspace)
        params = {
            'id':self.ws_id,
            'objects': [{
                'data': data,
                'name': objid,
                'type': "KBasePhenotypes.PhenotypeSet",
                'meta': {},
                'provenance': self.provenance()
            }]
        }
        self.ws_client().save_objects(params)
        self.obj_created.append({"ref":self.create_ref(objid,self.ws_name),"description":""})

    def save_solution_as_fba(self,fba_or_solution,mdlutl,media,fbaid,workspace=None,fbamodel_ref=None,other_solutions=None):
        if not isinstance(fba_or_solution,MSFBA):
            fba_or_solution = MSFBA(mdlutl,media,primary_solution=fba_or_solution)
        fba_or_solution.id = fbaid
        if other_solutions != None:
            for other_solution in other_solutions:
                fba_or_solution.add_secondary_solution(other_solution)
        data = fba_or_solution.generate_kbase_data(fbamodel_ref,media.info.reference)
        #If the workspace is None, then saving data to file
        if not workspace and self.util:
            self.util.save(fbaid,data)
        else:
            #Setting the workspace
            if workspace:
                self.set_ws(workspace)
            #Setting provenance and saving model using workspace API
            params = {
                'id':self.ws_id,
                'objects': [{
                    'data': data,
                    'name': fbaid,
                    'type': "KBaseFBA.FBA",
                    'meta': {},
                    'provenance': self.provenance()
                }]
            }
            self.ws_client().save_objects(params)
            self.obj_created.append({"ref":self.create_ref(fbaid,self.ws_name),"description":""})