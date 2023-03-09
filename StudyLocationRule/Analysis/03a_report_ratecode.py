import os
import logging
import numpy as np
import neuropype.nodes as nn
import neuropype.engine as ne
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'misc'))
from misc.misc import sess_infos
from functools import reduce


# Config
logging.basicConfig(level=logging.INFO)
preproc_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Data', 'Preprocessed'))
target_theta = np.deg2rad([270, 315, 0, 45, 90, 135, 180, 225])
target_str = np.array(['UU', 'UR', 'RR', 'DR', 'DD', 'DL', 'LL', 'UL'])  # The y-zero was up, so flip u/d.
target_mapping = {target_str[_]: theta for _, theta in enumerate(target_theta)}

slice_defs = [
    {'event': 'Target', 'range': (-0.60, 0.2)},  # stop time lengthened for longest target-cue onset asynchrony.
    {'event': 'Cue', 'range': (-0.05, 1.25)},
    {'event': 'Go', 'range': (-0.05, 0.30)}
]  # Total is Cue + Go (1.3 + 0.35 = 1.65) + longest target period
# which begins at 0.6 sec before target onset and ends at cue onset.

sac_rt_range = (-0.01, 0.60)  # ReactionTime is relative to Go cue. ReactionTime must be within this range.

TARGET_RATE = 100

# TODO: Why is sess_ix 3 returning
"""
Importing session sra3_1_j_052_00+
Found 1398 trials, 1652 timestamps(-0.6000 to 1.0510 at 1000.0 Hz), 34 units
"""


if __name__ == "__main__":
    # For testing, limit to 1 or a few sessions:
    # sess_infos = [_ for _ in sess_infos if _['exp_code'] in ['sra3_1_j_050_00+']]
    for sess_ix, sess_info in enumerate(sess_infos):
        if sess_ix < 0:
            # Skip   ^ files while debugging
            continue

        behav_pkt = nn.ImportH5(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_behav.h5'))()
        gaze_pkt = nn.ExtractStreams(stream_names=['gaze'])(data=behav_pkt)
        behav_pkt = nn.ExtractStreams(stream_names=['markers'])(data=behav_pkt)

        # Only keep markers/events from trials that were completed, either to correct target or incorrect target.
        occ_criteria = [{'name': 'newOutcomeCode', 'value': 0}, {'name': 'newOutcomeCode', 'value': 9},
                        {'name': 'OutcomeCode', 'value': 0}, {'name': 'OutcomeCode', 'value': 9}]
        behav_pkt = nn.SelectInstances(selection=occ_criteria, combine_selections='or',
                                       combine_previous='and')(data=behav_pkt)

        # Preview reaction times to help select a good range.
        mrk_iax = behav_pkt.chunks['markers'].block.axes['instance']
        b_go = mrk_iax.data['Marker'] == 'Go'
        rts = mrk_iax.data[b_go]['ReactionTime']
        print(f"ReactionTime range: {np.nanmin(rts)} - {np.nanmax(rts)}")
        if False:
            # Helper to visualize reaction time distribution.
            import matplotlib.pyplot as plt
            # plt.switch_backend('Qt5Agg')
            plt.hist(rts, 50)
            plt.show()
        selection_criteria = [{'name': 'ReactionTime', 'operator': 'greater_equal', 'value': sac_rt_range[0]},
                              {'name': 'ReactionTime', 'operator': 'less_equal', 'value': sac_rt_range[1]}]
        behav_pkt = nn.SelectInstances(selection=selection_criteria, combine_previous='and')(data=behav_pkt)

        # Sometimes we don't have enough trials of a given condition to even split it up to train/test for ML.
        # Drop all conditions that have fewer than 10 trials.
        idat = behav_pkt.chunks['markers'].block.axes['instance'].data
        idat = idat[idat['Marker'] == 'Go']
        class_ids, counts = np.unique(idat['TargetClass'], return_counts=True)
        sel_dicts = [{'name': 'TargetClass', 'value': tc} for ix, tc in enumerate(class_ids) if counts[ix] >= 10]
        if len(sel_dicts) > 0:
            behav_pkt = nn.SelectInstances(selection=sel_dicts, combine_selections='or')(data=behav_pkt)

        # Smooth and decimate the spike rates.
        spk_pkt = nn.ImportFile(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_spk.h5'))()
        dec_factor = int(np.ceil(spk_pkt.chunks['spikerates'].block.axes['time'].nominal_rate / TARGET_RATE))
        rates = nn.ExtractStreams(stream_names=['spikerates'])(data=spk_pkt)
        rates = nn.IIRFilter(frequencies=[40, 50], mode='lowpass', offline_filtfilt=True)(data=rates)
        rates = nn.Decimate(factor=dec_factor, axis='time')(data=rates)  # prev_rate --> ~100 Hz

        # Smooth and decimate gaze.
        gaze_pkt = nn.IIRFilter(frequencies=[40, 50], mode='lowpass', offline_filtfilt=True)(data=gaze_pkt)
        gaze_pkt = nn.Decimate(factor=10, axis='time')(data=gaze_pkt)  # 1 kHz --> 100 Hz

        # Load LFPs if they exist.
        try:
            lfp_fname = os.path.join(preproc_root, sess_info['exp_code'] + '_lfp.h5')
            lfp_pkt = nn.ImportFile(filename=lfp_fname)()
        except OSError:
            print("Did not find file for LFP packet.")
            import neuropype.engine as ne
            lfp_pkt = ne.Packet()

        # Join different data streams together and segment.
        full_pkt = nn.MergeStreams(replace_if_exists=True)(dataN=[lfp_pkt, spk_pkt, rates, behav_pkt, gaze_pkt])

        if False:
            # Load all the data into memory (gasp!).
            # TODO: The lazy slicing is still a bit problematic in ne.concat below.
            full_pkt.chunks['spikerates'].block._data = full_pkt.chunks['spikerates'].block.data
            full_pkt.chunks['gaze'].block._data = full_pkt.chunks['gaze'].block.data
            full_pkt.chunks['markers'].block._data = full_pkt.chunks['markers'].block.data

        if False and sess_ix == 1:
            # Manually slice out all the trials from block 8 with their full lengths
            # (must retain go-cue time separately)
            full_iax = full_pkt.chunks['markers'].block.axes['instance']
            b_full_targ = full_iax.data['Marker'] == 'Target'
            b_full_go = full_iax.data['Marker'] == 'Go'
            go_id_for_targ = np.searchsorted(np.where(b_full_go)[0], np.where(b_full_targ)[0])
            assert len(np.unique(go_id_for_targ)) == len(go_id_for_targ)
            full_iax_tvec = behav_pkt.chunks['markers'].block.axes['instance'].times
            t_targ = full_iax_tvec[b_full_targ]
            t_go = full_iax_tvec[b_full_go] - t_targ
            t_sacc = t_go + rts

            append_interval = np.zeros_like(full_iax.data['ReactionTime'])
            append_interval[b_full_targ] = t_go
            full_iax.append_fields(['TargGoInterval'], [append_interval], props=0)

            tr_dur = 2.0  # np.max(t_sacc)
            long_pkt = nn.Segmentation(time_bounds=(-0.6, tr_dur), select_markers=['Target'])(data=full_pkt)

            long_iax = long_pkt.chunks['spiketrains'].block.axes['instance']
            b_kept = np.in1d(full_iax.data['TrialIdx'], long_iax.data['TrialIdx'])
            assert np.array_equal(full_iax_tvec[np.logical_and(b_full_targ, b_kept)],
                                  long_pkt.chunks['spiketrains'].block.axes['instance'].times)

            long_fname = sess_info['exp_code'] + '_long.h5'
            long_fname = long_fname.replace("+", "")
            nn.ExportH5(filename=os.path.join(preproc_root, long_fname))(data=long_pkt)

        # Segmentation is a bit complicated because there is variability between target onset and cue onset,
        # but we are not interested in that variability. We do multiple slices and concatenate.
        # However, for the per-target slice we use a slice length that captures the slice in the longest trial,
        # and for the shorter trials (where this long slice would duplicate data in the cue period), we
        # overwrite the duplicated data with nan.
        mrk_iax = full_pkt.chunks['markers'].block.axes['instance']
        targ_idx = mrk_iax.data['Marker'] == 'Target'
        cue_idx = mrk_iax.data['Marker'] == 'Cue'
        targ_cue_oa = mrk_iax.times[cue_idx] - mrk_iax.times[targ_idx]
        max_oa = np.round(200 * np.max(targ_cue_oa)) / 200  # Rounded to nearest 5 msec.
        slice_defs[0]['range'] = (slice_defs[0]['range'][0], max_oa + slice_defs[1]['range'][0])
        seg_pkts = []
        for seg_ix, slice_def in enumerate(slice_defs):
            seg_pkts.append(nn.Segmentation(time_bounds=slice_def['range'],
                                            select_markers=[slice_def['event']])(data=full_pkt))

        min_common_trial_ids = reduce(
            np.intersect1d, [reduce(np.intersect1d,
                                    [_.chunks[chnk_name].block.axes['instance'].data['TrialIdx'] for _ in seg_pkts])
                             for chnk_name in ['spikerates', 'spiketrains', 'gaze']])

        # Segmentation might lose trials due to discontinuities in the data, and lost trials could vary between slices
        # So we drop trials that are not common to all segments.
        for seg_ix, seg_pkt in enumerate(seg_pkts):
            seg_pkts[seg_ix] = nn.SelectInstances(selection={'TrialIdx': min_common_trial_ids})(data=seg_pkt)

        # For each chunk, replace the samples at the end of seg 0 that are duplicated in seg 1 with nan.
        for chnk_name in ['spikerates', 'spiketrains', 'gaze']:
            t_targ = seg_pkts[0].chunks[chnk_name].block.axes['instance'].times[:, None]\
                     + seg_pkts[0].chunks[chnk_name].block.axes['time'].times[None, :]
            cue_0 = seg_pkts[1].chunks[chnk_name].block.axes['instance'].times\
                    + seg_pkts[1].chunks[chnk_name].block.axes['time'].times[0]
            b_dup = t_targ > cue_0[:, None]
            b_dup = np.repeat(b_dup[:, :, None], seg_pkts[0].chunks[chnk_name].block.shape[-1], axis=2)
            seg_pkts[0].chunks[chnk_name].block.data[b_dup] = np.nan

        # Now we need to adjust the timestamps of the not-first segments so that they follow each other.
        # This has to be done chunk-by-chunk because of the potentially different sampling rates.
        out_pkt = seg_pkts[0]
        for seg_ix in range(1, len(seg_pkts)):
            for chnk_name in seg_pkt.chunks:
                out_blk = out_pkt.chunks[chnk_name].block
                seg_blk = seg_pkts[seg_ix].chunks[chnk_name].block
                # Fix the time vector
                new_t0 = out_blk.axes['time'].times[-1] + (1 / out_blk.axes['time'].nominal_rate)
                old_t0 = seg_blk.axes['time'].times[0]
                seg_blk.axes['time'].times = seg_blk.axes['time'].times - old_t0 + new_t0
                # Concatenate. All instance data should be the same so this SHOULD work.
                out_pkt.chunks[chnk_name].block = ne.concat(ne.time, out_blk, seg_blk)
            seg_pkts[seg_ix] = None  # Save memory

        # Save the segmented data.
        # Kaggle strips out the "+" symbol from the filename so we'll do that preemptively.
        seg_fname = sess_info['exp_code'] + '_segmented.h5'
        seg_fname = seg_fname.replace("+", "")
        nn.ExportH5(filename=os.path.join(preproc_root, seg_fname))(data=out_pkt)

        if False:
            psth_edges = [(-0.25, 0), (0, 0.25), (0.25, 1.25), (1.25, np.inf)]
            psth_pkt = nn.SelectRange(axis='instance-fields', selection=['TargetStr'], unit='names')(data=seg_pkt)
            psth_pkt = nn.GroupedMean(error_type='sem')(data=psth_pkt)
            psth_pkt = nn.AssignTargets(mapping=target_mapping, iv_column='TargetStr')(data=psth_pkt)
            rates_only = nn.ExtractStreams(stream_names=['spikerates'])(data=psth_pkt)
            nn.PSTHRadarPlot(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_psth.pdf'), units_per_page=4,
                             window_edges=psth_edges)(data=rates_only)
