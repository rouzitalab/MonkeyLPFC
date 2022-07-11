"""
author: chadwick.boulay@gmail.com

"""

import os
import logging
import numpy as np
from neuropype.engine import Graph, CPE
import neuropype.nodes as nn

from custom_nodes.ReshapeDPCA import ReshapeDPCA


logging.basicConfig(level=logging.INFO)

# Define data pairs.
data_root = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
session_names = ['sra3_1_j_050_00+']

# TODO: Optimal Linear Estimator: http://www.neuro.uni-bremen.de/~dip/estimator.html
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2783655/

for sess_name in session_names:
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '.h5')),
                nn.ExtractStreams(stream_names=["spikerates", "segmented-markers"]))
    cpe = CPE(graph=patch)
    cpe.loop_run()
