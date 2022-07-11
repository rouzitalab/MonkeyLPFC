"""
author: chadwick.boulay@gmail.com

"""

import os
import logging
import numpy as np
from neuropype.engine import Graph, CPE
import neuropype.nodes as nn


logging.basicConfig(level=logging.INFO)

# Define data pairs.
data_root = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
session_names = ['sra3_1_j_050_00+']

for sess_name in session_names:
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '.h5')),
                nn.ExtractStreams(stream_names=["spikerates", "segmented-markers"]),
                nn.Decimate(factor=20, axis='time'),  # 1 kHz --> 50 Hz
                nn.GroupedMean(),
                nn.ExportH5(filename=os.path.join(data_root, sess_name, sess_name + '_psth.h5')))
    cpe = CPE(graph=patch)
    cpe.loop_run()

for sess_name in session_names:
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '_psth.h5')),
                nn.PSTHRadarPlot(filename=sess_name + '_psth.pdf', units_per_page=4,
                                 output_root=os.path.join(data_root, sess_name),
                                 window_edges=[(-np.inf, 0), (0, 0.5), (0.5, 1.5), (1.5, np.inf)]))
    cpe = CPE(graph=patch)
    cpe.loop_run()
