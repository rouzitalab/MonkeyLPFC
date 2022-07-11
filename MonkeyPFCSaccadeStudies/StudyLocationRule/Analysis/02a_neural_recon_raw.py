"""
author: chadwick.boulay@gmail.com

This script processes neural data obtained from 32-channel arrays in macaque PFC while
they performed a delayed saccade task.
The neural data are in Blackrock NSX and NEV format and are imported with the ImportNSX node.

The script aims to give a reconstruction of high-frequency raw data for the purposes of offline spike sorting.
If we do not have 30 kHz full-band data:
    Spike waveforms are superimposed on zeros at the original 30 kHz.
    If available, LFPs are upsampled to 30 kHz and superimposed on the spike waveforms.


"""

import os
import logging
import numpy as np
import json
import neuropype.nodes as nn
import neuropype.engine as ne
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'misc'))
import misc
from misc.misc import sess_infos, map_files
import custom_neuropype as cn


# Config
logging.basicConfig(level=logging.INFO)
data_root = misc.get_data_root()
preproc_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Data', 'Preprocessed'))

ADD_NOISE = False
ADD_LFPS = False


if __name__ == "__main__":
    for sess_ix, sess_info in enumerate(sess_infos):
        if sess_ix < 0:
            # Skip   ^ files while debugging
            continue

        # Load the datafile and do minimal processing
        dat_filename = os.path.join(data_root, 'Cerebus_DATA',
                                    sess_info['name'], sess_info['date'], sess_info['nsx'])
        nsx_pkt = nn.ImportNSX(filename=dat_filename, load_signals=False)()  # 'events', 'spiketrains', 'waveforms'

        if 'analogsignals' in nsx_pkt.chunks:  # <-- False if ^^ load_signals=False
            # Rename elec{N} to ch{N}
            nsx_pkt = nn.RenameChannels(streams=['analogsignals'], search_expr='elec', replace_expr='ch')(data=nsx_pkt)

            # Drop the Reward channel. analog_pkt can then be converted to int16
            analog_pkt = nn.ExtractStreams(stream_names=['analogsignals'])(data=nsx_pkt)
            ch_strs = ['ch' + str(_) for _ in range(256)]
            keep_chans = [_ for _ in analog_pkt.first().block.axes['space'].names if _ in ch_strs]
            analog_pkt = nn.SelectRange(axis='space', selection=keep_chans, unit='names')(data=analog_pkt)
            analog_pkt.chunks['analogsignals'].block.data[:] =\
                analog_pkt.chunks['analogsignals'].block.data.astype(np.int16)
            nsx_pkt = nn.MergeStreams(replace_if_exists=True)(dataN=[nsx_pkt, analog_pkt])
            del ch_strs

        nsx_pkt = misc.ParseDigitalWords()(data=nsx_pkt)

        # Get electrode positions from .cmp file.
        map_file = os.path.join(data_root, 'Cerebus_DATA', sess_info['name'], map_files[sess_info['name']])
        nsx_pkt = cn.UpdateElectrodePositions(filename=map_file, banks=[sess_info['bank']])(nsx_pkt)

        # Generate continuous signals from waveforms. Note that this collapses sorting.
        # TODO: Use DatasetView and fix slicing
        recon_pkt = cn.ReconContFromEvents(add_noise=False, add_lfps=False)(data=nsx_pkt)

        nsx_pkt = nn.ExtractStreams(stream_names=['spiketrains', 'events'])(data=nsx_pkt)
        nsx_pkt = nn.MergeStreams(replace_if_exists=True)(dataN=[nsx_pkt, recon_pkt])

        # Specify chunk size of 4 channels x 25 MB, and chunk cache of 4 GB
        chunk_data = {
            'events': True,
            'spiketrains': True,
            'recon_raw': (4, 25*1024*1024)
        }
        nn.ExportH5(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_recon_raw.h5'),
                    chunk_data=chunk_data, compression='gzip',
                    rdcc_nbytes=4000*1024*1024)(data=nsx_pkt)
