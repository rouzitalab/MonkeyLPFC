"""
author: chadwick.boulay@gmail.com

"""

import os
import logging
from neuropype.engine import Graph, CPE
import neuropype.nodes as nn


logging.basicConfig(level=logging.INFO)

# Define data pairs.
data_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                         '..', 'data'))
session_names = ['sra3_1_j_050_00+']

# Spike-LFP coherence pype
for sess_name in session_names:
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '.h5')),
                nn.ExtractStreams(stream_names=["segmented-markers", "spiketrains"]),
                nn.SegmentedSignalPlot())
    cpe = CPE(graph=patch)
    cpe.loop_run()
