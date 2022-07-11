"""
author: chadwick.boulay@gmail.com

This script processes neural data obtained from 32-channel arrays in macaque PFC while
they performed a delayed saccade task.
The neural data are in Blackrock NSX and NEV format and are imported with the ImportNSX node.
This yields 4 chunks: analogsignals (1 kHz), events, spiketrains, waveforms.

Events are processed from the Blackrock data and compared with the previously processed behavioural data
(see 01_preprocess_behav.py). Event alignment is used to calcluate the timestamp offset between the data types
and the neural data clock is changed to be the same as the behavioural data clock.

We then do processing that can be calculated from unsegmented data.
1. Clean spike trains. a - undo sorting, b - sanitize, c - drop waveforms

Subsequent analysis scripts can load this hdf5 file, use (AssignTargets and) Segmentation,
then calculate trial-based stats.

then use an AssignTargets node to select one marker type (e.g. Cue-) as the time-locking event,
and to assign the value of the target location to the event.
Once the markers are set, it is then possible to use Segmentation node to segment the data around the marker event.
After segmentation, it is possible to calculate per-trial features and visualize the results.

"""

import os
import logging
import neuropype.nodes as nn
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

EXPORT_LFP = False
DROP_NOISE_UNITS = True
DEBUG_SMALL = False
USE_CACHING = False          # Set to True when debugging. Otherwise False (and clear cache folder!)

if __name__ == "__main__":
    for sess_ix, sess_info in enumerate(sess_infos):
        if sess_ix < 1:
            # Skip   ^ files while debugging
            continue

        # Import the manual adjusted CortexPhy folder
        phy_pkt = nn.ImportCortexPhy(foldername=os.path.join(preproc_root, sess_info['exp_code'], 'phy'))()

        # Import the neural data
        sorted_fname = os.path.join(preproc_root, sess_info['exp_code'] + '_sorted.h5')
        nsx_pkt = nn.ImportH5(filename=sorted_fname)()  # events, recon_raw, waveforms
        # Separate the continuous time series data from nsx_pkt
        sig_pkt = nn.ExtractStreams(stream_names=['recon_raw'])(data=nsx_pkt)
        nsx_pkt = nn.ExtractStreams(stream_names=['events', 'waveforms'])(data=nsx_pkt)

        # Merge phy sorted info into nsx_packet.chunks['waveforms']
        phy_iax = phy_pkt.chunks['phy'].block.axes['instance']
        wf_iax = nsx_pkt.chunks['waveforms'].block.axes['instance']
        import numpy as np
        assert np.array_equal(wf_iax.data['time_ind'], phy_iax.data['time_ind'])
        wf_chan_id = np.array([int(_[2:]) for _ in wf_iax.data['chan_label']]) - 1
        assert np.array_equal(wf_chan_id,
                              phy_iax.data['chan_id'])
        # _, wf_iax_chan_ids = np.unique(wf_iax.data['chan_label'], return_inverse=True)
        # np.unique(1000*wf_iax_chan_ids + phy_iax.data['chan_id'])
        wf_iax.data['unit_id'] = phy_iax.data['unit_id']
        clust_prop = phy_iax.data.dtype.descr[list(phy_iax.data.dtype.names).index('clust_grp')][0][0]
        wf_iax.append_fields(['clust_grp'], field_data=[phy_iax.data['clust_grp']], props=[clust_prop])

        if DROP_NOISE_UNITS:
            from neuropype.engine import instance
            wf_iax = nsx_pkt.chunks['waveforms'].block.axes['instance']
            b_noise = wf_iax.data['clust_grp'] == 'noise'
            # Load the data into memory. Otherwise the integer slicing on the very long instance axis will be too slow.
            nsx_pkt.chunks['waveforms'].block._data = nsx_pkt.chunks['waveforms'].block.data
            nsx_pkt.chunks['waveforms'].block = nsx_pkt.chunks['waveforms'].block[..., instance[~b_noise], ...]

        # Get sparse event packet from waveform packet.
        spk_pkt = nn.CreateSparseEventsFromWaveforms(
            new_chunk_name='spiketrains', axes_source=sig_pkt, use_caching=USE_CACHING)(
            data=nn.ExtractStreams(stream_names=['waveforms'])(data=nsx_pkt))
        nsx_pkt = nn.MergeStreams(replace_if_exists=True)(dataN=[nsx_pkt, spk_pkt])
        spk_pkt = None
        sig_pkt = None  # recon_raw no longer needed.

        # Import the behaviour. It must have NEVTimestamps field already from 01_preprocess_behav.py
        behav_filename = os.path.join(preproc_root, sess_info['exp_code'] + '_behav.h5')
        behav_pkt = nn.ImportH5(filename=behav_filename)()

        # From the neural data, subtract the temporal offset between it and the behaviour timestamps.
        mrk_pkt = nn.ExtractStreams(stream_names=['markers'])(data=behav_pkt)
        tshifts = nn.GetTimeShifts(other_source='NEVTimestamps', return_inverse=False
                                   )(data=mrk_pkt, return_outputs='all')['offset']
        nsx_pkt = nn.ShiftTimestamps(offset=tshifts, include_markers=True)(data=nsx_pkt)
        behav_pkt = None
        mrk_pkt = None

        if DEBUG_SMALL:
            from neuropype.engine import time, instance
            # Make a smaller file for debugging.
            tvec = sig_pkt.chunks['recon_raw'].block.axes['time'].times
            tstop = tvec[0] + 180
            sig_pkt.chunks['recon_raw'].block = sig_pkt.chunks['recon_raw'].block[
                ..., time[tvec < tstop], ...]

            st_ts = nsx_pkt.chunks['spiketrains'].block.axes[time].times
            nsx_pkt.chunks['spiketrains'].block = nsx_pkt.chunks['spiketrains'].block[
                ..., time[st_ts < tstop], ...]

            wf_ts = nsx_pkt.chunks['waveforms'].block.axes[instance].times
            nsx_pkt.chunks['waveforms'].block = nsx_pkt.chunks['waveforms'].block[
                ..., instance[wf_ts < tstop], ...]

            ev_ts = nsx_pkt.chunks['events'].block.axes[instance].times
            nsx_pkt.chunks['events'].block = nsx_pkt.chunks['events'].block[
                ..., instance[ev_ts < tstop], ...]

        # Drop events and waveforms from nsx pkt before saving.
        if EXPORT_LFP and 'analogsignals' in nsx_pkt.chunks:
            lfp_pkt = nn.ExtractStreams(stream_names=['analogsignals'])(data=nsx_pkt)
            nn.ExportH5(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_lfp.h5'))(data=lfp_pkt)

        # Convert (a copy of) spiketrains/rasters to spikerates.
        spk_pkt = nn.ExtractStreams(stream_names=['spiketrains'])(data=nsx_pkt)
        spk_pkt = nn.SanitizeEventTrain(downsample_rate=1000)(data=spk_pkt)
        rate_pkt = nn.InstantaneousEventRate(kernel='alpha', kernel_parameter=0.05)(data=spk_pkt)
        rate_pkt = nn.RenameStreams(name_changes={'spiketrains': 'spikerates'})(data=rate_pkt)

        spk_pkt = nn.MergeStreams(replace_if_exists=True)(data1=spk_pkt, data2=rate_pkt)
        nn.ExportH5(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_spk.h5'))(data=spk_pkt)
