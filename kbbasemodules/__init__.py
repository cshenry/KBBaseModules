# -*- coding: utf-8 -*-

from __future__ import absolute_import

# set the warning format to be on a single line
import sys
import logging

__author__ = "Christopher Henry"
__email__ = "chenry@anl.gov"
__version__ = "0.0.1"

logger = logging.getLogger(__name__)

print("KBBaseModules", __version__)

import kbbasemodules
from kbbasemodules.clients import (
    DataFileUtil,
    RAST_SDK,
    KBaseReport,
    chenry_utility_module,
    AssemblyUtil,
    cb_annotation_ontology_api,
    GenomeFileUtil,
    AbstractHandle,
    Workspace
)