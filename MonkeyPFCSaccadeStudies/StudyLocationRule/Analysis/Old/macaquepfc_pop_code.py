"""
author: chadwick.boulay@gmail.com

"""

import os
import logging
import neuropype.nodes as nn
from neuropype.nodes.custom import DPCAPlot, CPDecompPlot


logging.basicConfig(level=logging.INFO)

# Define data pairs.
data_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                         '..', 'data'))
session_names = ['sra3_1_j_050_00+']

dpca_groups = [[0, 1, 2, 3],
               [4, 5, 6, 7]]  # * np.pi / 4
output_margs = ['target_pair', 'target_pair', 'decision']

# calculate and plot features
# Note that we do not use an intermediate h5 file here because much of the plotted information
# comes from the model stored in the metadata, which maybe isn't currently exporting correctly.
for sess_name in session_names:
    data_file = os.path.join(data_root, sess_name, sess_name + '_segs.h5')
    data_pkt = nn.ImportH5(filename=data_file)()
    data_pkt = nn.ExtractStreams(data=data_pkt, stream_names=["spikerates", "segmented-markers"])()
    data_pkt = nn.Decimate(data=data_pkt, factor=20, axis='time')()  # 1 kHz --> 50 Hz

    # dPCA
    dpca_pkt = nn.DemixingPCA(data=data_pkt,
                              grouping=dpca_groups,
                              labels='ds',  # first grouping dim is 'decision', second is 'stimulus'
                              join={'target_pair': ['s', 'st'], 'decision': ['d', 'dt'],
                                    'time': ['t'], 'interaction': ['ds', 'dst']},
                              component_sort=output_margs)()
    # TODO: Convert DPCAPlot to report. Use new 'model' output from DemixingPCA
    nn.DPCAPlot(data=dpca_pkt,
                filename=sess_name + '_dPCA2.pdf',
                output_root=os.path.join(data_root, sess_name),
                margs_3d=output_margs)()

    # CPD
    cpd_pkt = nn.TensorDecomposition(data=data_pkt,
                                     num_components=4, do_projection=False)()
    # TODO: Convert CPDecompPlot to report. Use new 'model' output.
    nn.TensorDecompositionPlot(data=cpd_pkt, filename=sess_name + '_cpdecomp2.pdf',
                               output_root=os.path.join(data_root, sess_name))()
