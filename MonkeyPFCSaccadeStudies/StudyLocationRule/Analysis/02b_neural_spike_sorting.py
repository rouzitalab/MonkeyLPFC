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
from neuropype.utilities.math import prev_prime
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'misc'))
import misc
from misc.misc import sess_infos
import custom_neuropype as cn


# Config
logging.basicConfig(level=logging.INFO)
data_root = misc.get_data_root()
preproc_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Data', 'Preprocessed'))

EXPORT_WF4AE = False        # Whether or not to export waveforms for autoencoder.
EXPORT_RECON_RAW = False    # Whether or not to export reconstructed raw timeseries
USE_CACHING = False         # Set to True when debugging. Otherwise False (and clear cache folder!)
DEBUG_SEG = False
RETHRESHOLD = True
WF_WIN = [-10 / 30e3, 38 / 30e3]


if __name__ == "__main__":
    for sess_ix, sess_info in enumerate(sess_infos):
        if sess_ix < 5:
            # Skip   ^ files while debugging
            continue

        # Load the (reconstructed) raw data, including events and spiketrains.
        recon_fname = os.path.join(preproc_root, sess_info['exp_code'] + '_recon_raw.h5')
        if DEBUG_SEG:
            seg_fname = os.path.join(preproc_root, sess_info['exp_code'] + '_recon_raw_seg.h5')
            if not os.path.exists(seg_fname):
                _pkt = nn.ImportH5(filename=recon_fname)()
                tvec = _pkt.chunks['recon_raw'].block.axes['time'].times
                seg_pkt = nn.SelectRange(selection=f"{tvec[0]+120}:{tvec[0]+180}",
                                         axis='time', unit='seconds')(data=_pkt)
                nn.ExportH5(filename=seg_fname)(data=seg_pkt)
            recon_fname = seg_fname

        # Because we have an expectation of what the dataset dimensions are, we can optimize
        # disk access a little by specifying per-dataset cache and number of slots
        rdcc_nbytes = 4000*1024*1024  # 4 GB. A 2 hr recording needs ~ 800 MB per channel, blocks are 4,X
        chunk_size = 100*1024*1024  # 100 MB - known from 02a.
        rdcc_nslots = 100 * max(1, int(rdcc_nbytes / chunk_size))
        rdcc_nslots = prev_prime(rdcc_nslots)
        _pkt = nn.ImportH5(filename=recon_fname, rdcc_nbytes=rdcc_nbytes, rdcc_nslots=rdcc_nslots)()

        sig_pkt = nn.ExtractStreams(stream_names=['recon_raw'])(data=_pkt)
        _pkt = nn.ExtractStreams(stream_names=['events', 'spiketrains'])(data=_pkt)

        # Convert data to float32 - necessary for whitening.
        sig_pkt = cn.AsType(dtype='float32', use_caching=USE_CACHING)(data=sig_pkt)

        # NODO: Highpass filter - only needed if added white noise and/or LFPs during reconstruction.

        # Drop samples where all channels have value==0
        trim_pkt = cn.GetNonzeroData(waveform_window=WF_WIN, use_caching=USE_CACHING)(
            data=nn.MergeStreams()(
                dataN=[sig_pkt, nn.ExtractStreams(stream_names=['spiketrains'])(data=_pkt)]))

        # Calculate spatial whitening matrices
        trim_pkt = nn.SignalWhitening(init_on=[0, 5.0], emit_calib_data=False,
                                      store_matrices=True, offline_use_all=True,
                                      offline_cores=1, regularization=1.0e-06,
                                      segsize=-1,
                                      use_caching=USE_CACHING)(data=trim_pkt)
        # TODO: Maybe whitening_mats should come as Blocks
        import numpy as np
        from neuropype.engine import Packet, Block, Chunk
        whitening_mats = trim_pkt.chunks['recon_raw'].props['whitening']  # sph_matrix and cov_matrix
        sph_matrix = np.array(whitening_mats['sph_matrix'])
        wht_pkt = Packet(chunks={'sphering': Chunk(
            block=Block(data=sph_matrix,
                        axes=(trim_pkt.chunks['recon_raw'].block.axes['space'],
                              trim_pkt.chunks['recon_raw'].block.axes['space']))
        )})
        trim_dat = None
        trim_pkt = None

        # Apply whitening to full signal
        # data chunking should be 4 channels x 26214400 samples (100 MB).
        # If we get a single time-chunk for all 32 channels, that's 838860800 samples (3.2 GB @ 32-bit).
        # We'll make our segsize an even 1/N fraction of that. If correct, we should see the first
        # fetch go slowly, then the next < N fetches go quickly because the data are already in memory.
        sig_pkt = nn.ApplyLinearTransformation(segsize=26214400 * 32 // 8,
                                               use_caching=USE_CACHING)(data=sig_pkt, filters=wht_pkt)

        # If SignalWhitening is not used then RectangularMovingWindowStandardization
        # or RobustRectangularStandardization is required for thresholds in next node to work properly.

        # Detect threshold crossing events
        if RETHRESHOLD:
            evt_pkt = nn.DetectSparseEvents(thresholds=[-5.0], align_at='peak',
                                            minimum_event_duration=60e-6,
                                            preserve_amplitude=True,
                                            new_chunk_name='spiketrains',
                                            use_caching=USE_CACHING)(data=sig_pkt)  # , adjacency_radius=0.0005
            print(f"Detected {len(evt_pkt.chunks['spiketrains'].block._data.indices)}. Original had "
                  f"{len(_pkt.chunks['spiketrains'].block._data.indices)}")
        else:
            evt_pkt = nn.ExtractStreams(stream_names=['spiketrains'])(data=_pkt)  # Note: Not peak-aligned.

        # Clean event trains: Drop the smaller of overlapping events
        evt_pkt = nn.SanitizeEventTrain(overlap_seek_window=300e-6,
                                        chan_pcnt_noise_thresh=30,
                                        offline_min_event_rate=0.1,
                                        min_window=WF_WIN,
                                        use_caching=USE_CACHING)(data=evt_pkt)

        # Get waveforms
        sig_pkt = nn.MergeStreams(replace_if_exists=True)(dataN=[sig_pkt, evt_pkt])
        wf_pkt = nn.ExtractEventWaveforms(time_bounds=WF_WIN,
                                          use_caching=USE_CACHING)(data=sig_pkt)

        if not RETHRESHOLD:
            wf_pkt = nn.AlignWaveformPeaks(align_upsample_factor=2,
                                           search_window=[-3 / 30e3, 3 / 30e3],  # +/- 3 samps
                                           use_caching=USE_CACHING)(data=wf_pkt)

        # Prepare for sorting - Calculate features
        # Some of the feature vectors might be nan if the waveform extended beyond data window,
        # though this should be handled by SanitizeEventTrain.
        ft_pkt = nn.GetWaveformFeatures(feature_types=["pca"],
                                        n_components=3,
                                        feature_fraction=1.0,
                                        use_caching=USE_CACHING)(data=wf_pkt)

        # Sorting via Clustering
        ft_pkt = nn.SimpleEventSorting(feature_fraction=1.0,
                                       branch_cluster_dimensions=3,
                                       cluster_method="isosplit5",
                                       branch_depth=2,
                                       use_caching=USE_CACHING)(data=ft_pkt)
        # Assign the sorted labels to the original waveforms
        wf_pkt.chunks['waveforms'].block.axes['instance'].data['unit_id'][:] = \
            ft_pkt.chunks['wf_feats'].block.axes['instance'].data['unit_id']

        # TODO: cluster_metrics
        # https://github.com/flatironinstitute/mountainsort/blob/master/packages/ms3/p_cluster_metrics.cpp
        # ? https://github.com/magland/ml_ephys/blob/master/ml_ephys/basic/p_compute_cluster_metrics.py
        #  and return isolation_metrics
        # https://github.com/flatironinstitute/mountainsort/blob/master/packages/ms3/p_isolation_metrics.cpp
        # TODO: ms3.combine_cluster_metrics
        # TODO: ms4alg.create_label_map
        # https://github.com/magland/ml_ms4alg/blob/master/ml_ms4alg/p_create_label_map.py
        # TODO: ms4alg.apply_label_map
        # https://github.com/magland/ml_ms4alg/blob/master/ml_ms4alg/p_apply_label_map.py

        print("Break here to plot some waveforms.")
        if False:
            import numpy as np
            import matplotlib.pyplot as plt

            wf_iax_dat = wf_pkt.chunks['waveforms'].block.axes['instance'].data
            test_chan = wf_iax_dat['chan_label'][np.argmax(wf_iax_dat['unit_id'])]
            test_chan = 'ch3'
            b_wfs = wf_iax_dat['chan_label'] == test_chan
            wf_dat = wf_pkt.chunks['waveforms'].block[b_wfs, :].data
            t_vec = wf_pkt.chunks['waveforms'].block.axes['time'].times * 1000
            labels = wf_iax_dat['unit_id'][b_wfs]
            uq_labels = np.unique(labels)

            from itertools import cycle
            colour_codes = map('C{}'.format, cycle(range(len(uq_labels))))
            class_colors = np.array([next(colour_codes) for _ in range(len(uq_labels))])

            for lab_ix, lab in enumerate(uq_labels):
                b_lab = labels == lab
                lab_dat = wf_dat[b_lab, :].T
                plt.plot(t_vec, lab_dat, color=class_colors[lab_ix], linewidth=0.2, zorder=0.0)
                plt.plot(t_vec, np.nanmean(wf_dat[b_lab], axis=0),
                         color=class_colors[lab_ix], linewidth=2.0, zorder=2.0)
            plt.show()

        # Merge results...
        _pkt = nn.MergeStreams(replace_if_exists=True)(dataN=[
            nn.ExtractStreams(stream_names=['events'])(data=_pkt),
            nn.ExtractStreams(stream_names=['recon_raw'])(data=sig_pkt),
            wf_pkt])  # sorted waveforms

        # Export to Phy
        from pathlib import Path
        parent = Path(preproc_root) / sess_info['exp_code']
        if not parent.is_dir():
            parent.mkdir()
        nn.ExportCortexPhy(output_path=os.path.join(preproc_root, sess_info['exp_code'], 'phy'))(
            data=_pkt, waveform_features=ft_pkt)

        # Export full rich dataset that cannot be reconstructed (easily) from phy
        nn.ExportH5(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_sorted.h5'),
                    chunk_data=True, compression='gzip')(data=_pkt)
