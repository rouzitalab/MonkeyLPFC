"""
author: chadwick.boulay@gmail.com

This script performs pre-processing and limited feature extraction
on data from 32-channel arrays in macaque PFC while they performed
a delayed saccade task.

The data are in Blackrock NSX format and are first imported with the ImportNSX node.
This yields 4 chunks: analogsignals (1 kHz), events, spiketrains, waveforms
Then the custom_nodes.FixEvents node is used to replace the 'events' chunk
with a 'markers' chunk that contains marker events for Target onset, Cue onset, and Saccade onset.
Note that the markers are loaded in from a .mat file accompanying each data file.

We then calculate features that make sense to calculate from unsegmented data.
1. Clean spike trains. a - undo sorting, b - sanitize, c - drop waveforms

We then save the data to an hdf5 file.

Subsequent analysis scripts can load this hdf5 file, use AssignTargets and Segmentation,
then calculate trial-based stats.


then use an AssignTargets node to select one marker type (e.g. Cue-) as the time-locking event,
and to assign the value of the target location to the event.
Once the markers are set, it is then possible to use Segmentation node to segment the data around the marker event.
After segmentation, it is possible to calculate per-trial features and visualize the results.

"""

import os
import logging
import numpy as np
from neuropype.engine import *
import neuropype.nodes as nn

from custom_nodes.FixEvents import FixEvents


logging.basicConfig(level=logging.INFO)


data_root = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                         '..', 'data'))
file_pairs = [('sra3_1_j_050_00+', 'datafile002.ns2'), ]

targ_angles = np.asarray([2, 1, 0, 7, 6, 5, 4, 3])  # * np.pi/4
targ_mapping = {'Target-' + str(class_id + 1): targ_angles[class_id] for class_id in range(8)}

for sess_name, data_file in file_pairs:
    data_pkt = nn.ImportNSX(filename=os.path.join(data_root, sess_name, data_file))()
    data_pkt = FixEvents(data=data_pkt,
                         filename=os.path.join(data_root, sess_name, sess_name + '.mat'))()
    data_pkt = nn.RenameChannels(streams=['analogsignals'], search_expr='elec', replace_expr='ch')(data=data_pkt)
    data_pkt = nn.AssignTargets(data=data_pkt, mapping=targ_mapping)()

    spk_pkt = nn.ExtractStreams(data=data_pkt, stream_names=['spiketrains', 'waveforms'])()
    spk_pkt = nn.SimpleSpikeSorting(data=spk_pkt, cluster_method="none")()
    spk_pkt = nn.ExtractStreams(data=spk_pkt, stream_names=['spiketrains'])()
    spk_pkt = nn.SanitizeSpikeTrain(data=spk_pkt, min_refractory_period=1.0,
                                    downsample_rate=1000, chan_pcnt_noise_thresh=30)()
    # spk_pkt = nn.InstantaneousEventRate(data=spk_pkt, kernel='gaussian', kernel_parameter=0.05)()
    # spk_pkt = nn.RenameStreams(data=spk_pkt, name_changes={'spiketrains': 'spikerates'})()

    data_pkt = nn.MergeStreams(data1=data_pkt, data2=spk_pkt, replace_if_exists=True)()
    data_pkt = nn.ExtractStreams(data=data_pkt,
                                 stream_names=['spiketrains', 'markers'])()
    nn.ExportH5(data=data_pkt, filename=os.path.join(data_root, sess_name,
                                                     sess_name + '_spiketrains.h5'))()
