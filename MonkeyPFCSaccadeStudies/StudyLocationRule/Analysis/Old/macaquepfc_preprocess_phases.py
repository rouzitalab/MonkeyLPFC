"""
Companion to macaquepfc_preprocess

"""

import os
import logging
import numpy as np

from neuropype.engine import Graph, CPE
import neuropype.nodes as nn

from custom_nodes.FixEvents import FixEvents



logging.basicConfig(level=logging.INFO)

data_root = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')
file_pairs = [('sra3_1_j_050_00+', 'datafile002.ns2'), ]

targ_angles = np.asarray([2, 1, 0, 7, 6, 5, 4, 3]) * np.pi/4
targ_mapping = {'Target-' + str(class_id + 1): targ_angles[class_id] for class_id in range(8)}
lfp_bands = [3, 8, 8, 18, 18, 30, 30, 55]

for sess_name, data_file in file_pairs:
    patch = Graph()
    patch.chain(nn.ImportNSX(filename=os.path.join(data_root, sess_name, data_file)),
                nn.FixEvents(filename=os.path.join(data_root, sess_name, sess_name + '.mat')),
                nn.RenameChannels(streams=['analogsignals'], search_expr='elec', replace_expr='ch'),
                nn.AssignTargets(mapping=targ_mapping),
                nn.ExtractStreams(stream_names=['analogsignals', 'markers']),
                nn.FilterBank(frequency_edges=lfp_bands, stop_atten=30, transition_widths=[2, ]),
                nn.ShiftTimestamps(compensate_filter_lag=True),
                nn.Segmentation(time_bounds=[-0.2, 1.799]),
                nn.Phase(),
                nn.ExtractStream(stream_name="analogsignals"),
                nn.RenameStreams(name_changes={"analogsignals": "phases"}),
                nn.ExportH5(filename=os.path.join(data_root, sess_name, sess_name + '_phases.h5')))

    cpe = CPE(graph=patch)
    cpe.loop_run()
