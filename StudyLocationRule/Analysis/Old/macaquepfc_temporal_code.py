"""
author: chadwick.boulay@gmail.com

"""

import os
import logging
from neuropype.engine import Graph, CPE
import neuropype.nodes as nn


logging.basicConfig(level=logging.INFO)

# Define data pairs.
data_root = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
session_names = ['sra3_1_j_050_00+']
lfp_bands = [(3, 8), (8, 18), (18, 30), (30, 55)]

# Spike-LFP coherence pype
for sess_name in session_names:
    importer = nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '.h5'))
    patch = Graph()
    patch.chain(importer,
                nn.ExtractStreams(stream_names=["analogsignals", "spiketrains", "segmented-markers"]),
                nn.MultitaperCoherence(half_bandwidth=5, ref_stream_name="analogsignals"),
                nn.ExtractStreams(stream_names=["spiketrains", "segmented-markers"]),
                nn.GroupedMean(),
                nn.ExportH5(filename=os.path.join(data_root, sess_name, sess_name + '_spklfpcoh.h5')))

# PAC pype
for sess_name in session_names:
    merger = nn.MergeStreams()
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '.h5')),
                nn.ExtractStreams(stream_names=["spikerates", "segmented-markers"]),
                (merger, 'data1'))
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '_phases.h5')),
                (merger, 'data2'))
    patch.chain(merger,
                nn.PhaseAmplitudeCoupling(phases_stream_name="phases", amplitudes_stream_name="spikerates"),
                nn.ExtractStreams(stream_names=["comodulogram", "segmented-markers"]),
                nn.GroupedMean(),
                nn.ExportH5(filename=os.path.join(data_root, sess_name, sess_name + '_spk_lfp_pac.h5')))
    cpe = CPE(graph=patch)
    cpe.loop_run()

for sess_name in session_names:
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '_spklfpcoh.h5')),
                nn.PSTHRadarPlot(filename=sess_name + '_spklfpcoh.pdf', units_per_page=4,
                                 output_root=os.path.join(data_root, sess_name),
                                 window_edges=lfp_bands))
    cpe = CPE(graph=patch)
    cpe.loop_run()

for sess_name in session_names:
    patch = Graph()
    patch.chain(nn.ImportH5(filename=os.path.join(data_root, sess_name, sess_name + '_spk_lfp_pac.h5')),
                nn.PSTHRadarPlot(filename=sess_name + '_spk_lfp_pac.pdf', units_per_page=4,
                                 output_root=os.path.join(data_root, sess_name),
                                 window_edges=lfp_bands))
    cpe = CPE(graph=patch)
    cpe.loop_run()
