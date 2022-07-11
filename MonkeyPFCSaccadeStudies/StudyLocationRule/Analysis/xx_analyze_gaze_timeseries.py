"""

"""

import os
import logging
import numpy as np
import platform
import neuropype.engine as ne
import neuropype.nodes as nn
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'misc'))
import misc


# Setup and constants
logging.basicConfig(level=logging.INFO)
TARG_ANGLES = np.asarray([2, 1, 0, 7, 6, 5, 4, 3])  # * np.pi/4
TARG_MAP = {'Target-' + str(_ + 1): TARG_ANGLES[_] for _ in range(8)}
data_roots = [
    os.path.abspath(r'E:\Chad\The Ottawa Hospital\Sachs Lab - Documents\RawData\Monkey'),  # SachsLab Tower - Windows
    os.path.abspath('/media/SachsLab/Chad/The Ottawa Hospital/Sachs Lab - Documents/RawData/Monkey'),  # SachsLab Tower - Linux
    os.path.abspath(r'D:\SachsLab\RawData\Common-Monkey'),  # Chad home - Windows
    os.path.abspath('/media/chad/STORE/SachsLab/RawData/Common-Monkey')  # Chad home - Linux
]
data_root = [_ for _ in data_roots if os.path.isdir(_)][0]
preproc_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Data', 'Preprocessed'))

# TODO: Iterate through files. For now we do one at a time.
sess_infos = [
    {'name': 'JerryLee', 'name_short': 'j', 'date': '090623', 'exp_code': 'sra3_1_j_050_00+', 'nsx': 'datafile002.ns2'}
]
sess_ix = 0
sess_info = sess_infos[sess_ix]

# Load the data from the psychophysics toolbox file (events)
data_pkt = nn.ImportFile(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_full.h5'))()

#  See Common_Matlab\Helper\getNewClass.m  --> Do it in analysis, not preprocessing

# LFPs
t = data_pkt.chunks['analogsignals'].block.axes[ne.time].times
srate = data_pkt.chunks['analogsignals'].block.axes[ne.time].nominal_rate
X = data_pkt.chunks['analogsignals'].block[:, ne.space[:-4]].data

# spiketrains
spiketrains = data_pkt.chunks['spiketrains'].block.data

# Gaze data
# TODO: Align 'gaze' with 'analogsignals' in time.
# y = data_pkt.chunks['gaze'].block[:, space[['PosX', 'PosY']]].data
# This y data is incorrect:
y = data_pkt.chunks['analogsignals'].block[:, ne.space[['xEye', 'yEye']]].data

#

import matplotlib.pyplot as plt
b_test = np.logical_and(t >= 3000.0, t < 3010.0)
plt.subplot(3, 1, 1)
plt.plot(t[b_test], X[b_test, :5])
plt.subplot(3, 1, 2)
plt.plot(t[b_test], y[b_test, :])
plt.subplot(3, 1, 3)
plt.eventplot(np.sort(t[spiketrains[0, b_test].nonzero()[1]]))
plt.show()

nn.InstantaneousEventRate
nn.FilterBank
nn.Magnitude
