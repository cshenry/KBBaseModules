from __future__ import absolute_import

import logging
import os
import sys
import json
import pandas as pd
from kbbasemodules.basemodule import BaseModule
from os.path import exists

logger = logging.getLogger(__name__)

class BaseAnnotationModule(BaseModule):
    def __init__(self,name,config,module_dir="/kb/module",working_dir=None,token=None,clients={},callback=None):
        BaseModule.__init__(self,name,config,module_dir,working_dir,token,clients,callback)
        self.genome_info_hash = {}
        self.annodata = {"Genomes":[],"Gene":[],"Type":[],"Function":[],"Score type":[],"Score":[]}
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
    
    def anno_client(self):
        if "cb_annotation_ontology_api" not in self.clients:
            from installed_clients.cb_annotation_ontology_apiClient import cb_annotation_ontology_api
            self.clients["cb_annotation_ontology_api"] = cb_annotation_ontology_api(self.callback_url,token=self.token)
        return self.clients["cb_annotation_ontology_api"]
    
    def genome_to_proteins(self,ref):
        output = self.get_object(ref,self.ws_id)
        self.genome_info_hash[ref] = output["info"]
        sequence_list = []
        for ftr in output["data"]["features"]:
            if "protein_translation" in ftr:
                sequence_list.append([ftr["id"],ftr["protein_translation"]])
        return sequence_list
    
    def process_genome_list(self,reference_list):
        genome_list = []
        for ref in reference_list:
            #TODO: add support for genomesets
            genome_list.append(ref)
        return genome_list
    
    def add_annotations_to_genome(self,genome_ref,suffix,annotations):
        """Loads specified gene annotation into KBase genome object
        
        Parameters
        ----------
        string - genome_ref
            KBase workspace reference to genome where annotations should be saved
        string - suffix
            Suffix to be used when saving modified genome back to KBase
        mapping<string gene_id,mapping<string ontology,mapping<string term,{"type":string,"score":float}>>> - annotations
            Annotations to be saved to genome
        Returns
        -------
        dict
            
        Raises
        ------
        """
        ontology_inputs = {}
        for geneid in annotations:
            for ontology in annotations[geneid]:
                if ontology not in ontology_inputs:
                    ontology_inputs[ontology] = {}
                if geneid not in ontology_inputs[ontology]:
                    ontology_inputs[ontology][geneid] = []
                for term in annotations[geneid][ontology]:
                    suffix = "" 
                    anno_data = {"term": term}
                    if "scores" in annotations[geneid][ontology][term]:
                        anno_data["metadata"] = [annotations[geneid][ontology][term]]
                        suffix = ";"+annotations[geneid][ontology][term]["scores"][1]+"="+str(annotations[geneid][ontology][term]["scores"][2])
                        self.annodata["Score"].append(annotations[geneid][ontology][term]["scores"][2])
                        self.annodata["Score type"].append(annotations[geneid][ontology][term]["scores"][1])
                    else:
                        self.annodata["Score"].append(None)
                        self.annodata["Score type"].append(None)
                    if "name" in annotations[geneid][ontology][term]:
                        anno_data["name"] = annotations[geneid][ontology][term]["name"]+suffix
                    else:
                        anno_data["suffix"] = suffix
                    self.annodata["Genomes"].append(genome_ref)
                    self.annodata["Gene"].append(geneid)
                    self.annodata["Type"].append(ontology)
                    self.annodata["Function"].append(term)
                    ontology_inputs[ontology][geneid].append(anno_data)
                        
        anno_api_input = {
            "input_ref":genome_ref,
            "output_name":self.get_genome_info(genome_ref)[1]+suffix,
            "output_workspace":self.ws_id,
            "overwrite_matching":1,
            "save":1,
            "provenance":self.provenance(),
            "events":[]
        }
        for ontology in ontology_inputs.keys():
            anno_api_input["events"].append({
                "ontology_id":ontology,
                "method":self.name+"."+self.method,
                "method_version":self.config["version"],
                "timestamp":self.timestamp,
                "ontology_terms":ontology_inputs[ontology]
            })
        anno_api_output = self.anno_client().add_annotation_ontology_events(anno_api_input)
        self.obj_created.append({"ref":anno_api_output["output_ref"],"description":"Saving PDB annotation for "+self.genome_info_hash[ref][1]})
        return anno_api_output