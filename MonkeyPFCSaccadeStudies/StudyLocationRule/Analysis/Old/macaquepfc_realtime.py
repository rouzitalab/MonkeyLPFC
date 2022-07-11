"""
author: chadwick.boulay@gmail.com


"""

import os
import logging
from qtpy import QtGui, QtCore, QtWidgets
import numpy as np
from neuropype.engine import *
import neuropype.nodes as nn


logging.basicConfig(level=logging.INFO)

data_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                         '..', 'data'))
file_pairs = [('sra3_1_j_050_00+', 'datafile002.ns2'), ]
targ_angles = np.asarray([2, 1, 0, 7, 6, 5, 4, 3])  # * np.pi/4
dpca_groups = [[0, 1, 2, 3],
               [4, 5, 6, 7]]  # * np.pi / 4
grp_strs = [str(x * np.pi/4) for x in targ_angles]

sess_ix = 0
sess_name = file_pairs[sess_ix][0]
data_file = file_pairs[sess_ix][1]

# file names
calib_file = os.path.join(data_root, sess_name, sess_name + '_segs.h5')
stream_file = os.path.join(data_root, sess_name, sess_name + '_spiketrains.h5')

# merge and split nodes need to be defined a priori
stream_node = nn.StreamData()
inject_node = nn.InjectCalibrationData()  # merge point
logreg_node = nn.LogisticRegression()  # split point

patch = Graph()

# calib data is already preprocessed. Just import and run into inject node.
patch.chain(
    nn.ImportH5(filename=calib_file),
    nn.ExtractStreams(stream_names=['spikerates', 'segmented-markers']),
    (inject_node, 'calib_data')
)

# Test streaming data needs a little processing first before it joints inject_node.
patch.chain(
    nn.ImportH5(filename=stream_file),
    stream_node,
    nn.RenameStreams(name_changes={'spiketrains': 'spikerates'}),
    nn.InstantaneousEventRate(kernel='gaussian', kernel_parameter=0.05),
    nn.Segmentation(time_bounds=[-0.2, 1.799], online_epoching='sliding'),
    (inject_node, 'streaming_data')
)

patch.chain(
    stream_node,
    nn.SpikeTrainRasterPlot()
)
# Extract features, classify
patch.chain(
    inject_node,
    nn.Decimate(factor=20, axis='time'),  # 1 kHz --> 50 Hz
    # DemixingPCA(grouping=dpca_groups,
    #             labels='ds',  # first grouping dim is 'decision', second is 'stimulus'
    #             n_components={'all': 5},
    #             join={'target_pair': ['s', 'st'], 'decision': ['d', 'dt'],
    #                   'time': ['t'], 'interaction': ['ds', 'dst']},
    #             component_sort=['target_pair', 'target_pair', 'decision']),
    # SelectRange(axis='time', selection='-40:'),
    nn.Mean(axis='time'),
    logreg_node,
    nn.MeasureLoss(accumulate_offline=True)
)

patch.chain(
    logreg_node,
    nn.OverrideAxis(old_axis='feature', new_axis='feature', init_data=grp_strs),
    # OverrideFeatureUnits(features_to_override=grp_strs, unit='radians', accept_space_axis=False),
    nn.BarPlotPG(x_tick_type='radians', x_label='Target', y_label='P(target)', bar_width=0.8, autoscale=False)
)

cpe = CPE(graph=patch)
# cpe.loop_run()


def update():
    global cpe
    cpe.loop_step()


if __name__ == '__main__':
    import sys
    app = QtWidgets.QApplication([])
    timer = QtCore.QTimer()
    timer.timeout.connect(update)
    timer.start(0)
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtWidgets.QApplication.instance().exec_()
