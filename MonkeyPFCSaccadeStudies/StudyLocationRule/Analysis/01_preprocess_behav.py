"""
author: chadwick.boulay@gmail.com

This script performs pre-processing on behavioural data obtained
from macaques while they performed a delayed saccaded task.

The event data are loaded from the psychophysics toolbox file with the ImportPTB node.
This creates a 'markers' chunk.

The raw eyetracker data are imported from the eyelink EDF file with the ImportEyelink node.
Note that this node is a work-in-progress. It needs better cross-platform unification.
This creates 'gaze' and 'gaze_events' chunks.

Note that the gaze_calib_adjust metadata from the ptb 'markers' chunk can be passed
into ImportEyelink and the data will be adjusted for calibration updates as it is imported.

We then save the data to an hdf5 file.

"""

import os
import logging
import numpy as np
import neuropype.nodes as nn
import neuropype.engine as ne
import matplotlib
import matplotlib.pyplot as plt
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..', 'misc'))
import misc
import custom_neuropype as cn
from misc.misc import sess_infos


matplotlib.interactive(True)
plt.switch_backend('Qt5Agg')
logging.basicConfig(level=logging.INFO)


# Dataset to load and preprocess
data_root = misc.get_data_root()
preproc_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'Data', 'Preprocessed'))


def debug_plot_gaze(pkt, trial_inds=None, n_rand_trials=5, plot_window=(-0.1, 0.4),
                    marker_ts_use_eyelink=False, timelock_event='Go'):
    """
    Plot a couple segments of eye tracking data
    :param pkt:
    :param trial_inds:
    :param n_rand_trials: Choose up to this many random trials if trial_inds is not provided.
    :param plot_window:
    :param marker_ts_use_eyelink: Set to True to use marker times from instance axis EyelinkTime, else it will use
        the regular marker timestamps. If using 'eyelink', then the continuous gaze timestamps (in 'gaze' block time
        axis) and the optional 'gaze_events' instance axis times should both use the eyelink clock. Otherwise, these
        should both use the same clock as marker timestamp clock (i.e. the PTB clock).
    :param timelock_event: one of 'Start', 'Fixation', 'Target', 'Cue', 'Delay', 'Go', 'Reward', 'Stop'
    :return:
    """

    if 'gaze' not in pkt.chunks:
        raise KeyError("gaze chunk required (until PTB packet includes gaze data?)")

    m_iax = pkt.chunks['markers'].block.axes[ne.instance]
    b_ev = m_iax.data['Marker'] == timelock_event

    if trial_inds is None:
        import random
        trial_inds = random.sample(m_iax.data[b_ev]['TrialIdx'].tolist(), n_rand_trials)
    trial_inds.sort()

    b_inds = np.in1d(m_iax.data['TrialIdx'], trial_inds)
    m_iax = pkt.chunks['markers'].block[..., ne.instance[b_ev & b_inds], ...].axes[ne.instance]
    m_times = m_iax.data['EyelinkTime'] if marker_ts_use_eyelink else m_iax.times

    if 'gaze_events' in pkt.chunks:
        e_blk = pkt.chunks['gaze_events'].block
        b_sac = e_blk.axes[ne.instance].data['Marker'] == 'Saccade'
        sac_times = e_blk[b_sac].axes[ne.instance].times
    else:
        sac_times = []

    g_blk = pkt.chunks['gaze'].block
    t_vec = np.arange(plot_window[0], plot_window[1], 1/g_blk.axes[ne.time].nominal_rate)  # For plotting
    t_inds = np.arange(len(t_vec)) - np.argmin(np.abs(t_vec))                              # For slicing
    for ix, idx in enumerate(trial_inds):
        t_0_ix = np.nanargmin(np.abs(g_blk.axes['time'].times - m_times[ix]))
        t_slice = t_inds + t_0_ix

        plt.subplot(len(trial_inds), 1, 1 + ix)
        for ch_ix, ch_name in enumerate(['gx_left', 'gy_left']):
            plt.plot(t_vec, g_blk[ne.time[t_slice], ne.space[ch_name]].data, label=ch_name)
            plt.ylabel("Trial {}".format(idx))
            plt.axvline(0, color='k', linestyle='--')
        if ix == len(trial_inds) - 1:
            plt.legend()
            plt.xlabel('Time after ' + timelock_event)

        b_trial_sac = np.logical_and(sac_times >= (m_times[ix] + t_vec[0]),
                                     sac_times <= (m_times[ix] + t_vec[-1]))
        if np.any(b_trial_sac):
            t_sac_times = sac_times[b_trial_sac] - m_times[ix]
            for t_sac_time in t_sac_times:
                plt.axvline(t_sac_time, color='r', linestyle=':')

        if 'newOutcomeCode' in m_iax.data.dtype.names:
            noc = m_iax.data[ix]['newOutcomeCode']
            txt = 'Target' if noc == 0 else ('Distractor' if noc == 9 else 'Unknown')
            plt.annotate(txt, (0, 0))
            plt.axhline(m_iax.data[ix]['sacEndX'], color='b', linestyle='--')
            plt.axhline(m_iax.data[ix]['sacEndY'], color='r', linestyle='--')

    plt.show()


if __name__ == "__main__":
    for sess_ix, sess_info in enumerate(sess_infos):
        if sess_ix < 0:
            # Skip   ^ files while debugging
            continue

        # Load the psychophysics toolbox file and (slowly) process the events.
        # Gives 'markers' chunk.
        ptb_filename = os.path.join(data_root, 'EDF and PTB files',
                                    sess_info['name_short'], 'ptbmat', sess_info['exp_code'] + '.ptbmat')

        mrk_pkt = misc.ImportPTB(filename=ptb_filename, use_caching=True)()

        # Load the Eyelink data
        eye_filename = os.path.join(data_root, 'EDF and PTB files',
                                    sess_info['name_short'], 'edf', sess_info['exp_code'] + '.edf')
        calib_updates = nn.ExtractChunkProperty(chunk_name_prefix='markers',
                                                property_name='gaze_calib_adjust')(data=mrk_pkt)
        eye_pkt = nn.ImportEyelink(filename=eye_filename, calib_updates=calib_updates)()

        # Convert gaze data units from pixels to degrees.
        ptb_params = nn.ExtractChunkProperty(chunk_name_prefix='markers', property_name='ptb_params')(data=mrk_pkt)
        eye_pkt = nn.SetGazeMetaData(screen_distance=ptb_params['subjectScreenDistance'] / 1000,
                                     screen_resolution=list(ptb_params['subjectScreenResolution']),
                                     screen_dimensions=list(ptb_params['subjectScreenSize'] / 1000))(data=eye_pkt)
        eye_pkt = nn.SelectRange(axis='space', selection=['gx_left', 'gy_left'], unit='names')(data=eye_pkt)
        eye_pkt = nn.GazeScreenCoordinateConversion(convert_from='pixels', convert_to='degrees',
                                                    convert_instance_fields=['gavx', 'gavy',
                                                                             'genx', 'geny',
                                                                             'gstx', 'gsty'])(data=eye_pkt)

        # full_pkt = nn.MergeStreams(replace_if_exists=True)(data1=mrk_pkt, data2=eye_pkt)
        # debug_plot_gaze(full_pkt, trial_inds=[1235,  218, 1401, 1104,  536], marker_ts_use_eyelink=True)

        # Shift EyeLink timestamps so they are on the same clock as the ptb timestamps. e.g., from ~2000 to about 17000.
        # EyeLink timestamps are discontinuous so we should shift each timestamp independently.
        tshifts = nn.GetTimeShifts(other_source='EyelinkTime')(data=mrk_pkt, return_outputs='all')
        eye_pkt = nn.ShiftTimestamps(timestamp_pairs=tshifts['timestamp_pairs'],
                                     include_markers=True, marker_fields=['end'])(data=eye_pkt)

        # Drop gaze samples with nan timestamps (outside of interpolation range)
        eye_pkt = nn.SelectRange(axis='time', selection='0...9999999999999', unit='seconds')(data=eye_pkt)

        # Replace blinks, _invalid, and _offscreen with nans. EyeLink has blink events.
        # Also, PTB? has [DEF.EYE_State_blink _invalid _offscreen]
        gaze_blk = eye_pkt.chunks['gaze'].block
        gaze_iax = eye_pkt.chunks['gaze_events'].block.axes['instance']
        blink_inds = np.where(gaze_iax.data['blink'])[0]
        for b_ix in blink_inds:
            b_blink = np.logical_and(gaze_blk.axes['time'].times > gaze_iax.times[b_ix],
                                     gaze_blk.axes['time'].times < gaze_iax.data['end'][b_ix])
            eye_pkt.chunks['gaze'].block.data[b_blink, :] = np.nan
        # Off-screen? Undetected blinks? Still have huge values after blinks.
        eye_pkt.chunks['gaze'].block.data[np.any(gaze_blk.data > 500, axis=1), :] = np.nan

        eye_pkt = nn.SanitizeGazeData()(data=eye_pkt)

        # Saccade/Fixation detector.
        gaze_events = nn.ExtractStreams(stream_names=['gaze'])(data=eye_pkt)
        gaze_events = cn.NSLRHMM()(data=gaze_events)
        gaze_events = nn.RenameStreams(name_changes={'gaze': 'gaze_events'})(data=gaze_events)

        # debug_plot_gaze(nn.MergeStreams(replace_if_exists=True)(data1=mrk_pkt, data2=eye_pkt, data3=gaze_events),
        #               # trial_inds=[1880, 1881, 1900, 1905],
        #                 plot_window=(-1.25, 0.4))

        # Merge saccade events into marker table
        sac_evs = nn.SelectInstances(selection=[{'name': 'Marker', 'value': 'Saccade'}])(data=gaze_events)
        mrk_pkt = nn.JoinIVTables(handle_common='assert-equal')(data=mrk_pkt, data2=sac_evs)
        mrk_pkt = nn.FillInstanceDataTable(fill_direction='ffill', fill_fields=['TrialIdx'])(data=mrk_pkt)

        # Calculate reaction time and save as new field. Make sure all events in trial have same value for ReactionTime.
        mrk_pkt = nn.CalculateMarkerTimeDifference(source_event='Go', target_event=['Saccade', 'Gaze'],
                                                   new_field_name='ReactionTime', look_for='nearest',
                                                   grouping_fields=['TrialIdx'])(data=mrk_pkt)
        mrk_pkt = nn.FillInstanceDataTable(fill_direction='ffill', grouping_fields=['TrialIdx'],
                                           fill_fields=['ReactionTime'])(data=mrk_pkt)
        # Fill saccade-preceding events in trial with saccade-related variables.
        trial_fill_fields = ['ReactionTime', 'StartX', 'StartY', 'TargetX', 'TargetY', 'DistractorX', 'DistractorY',
                             'PosX', 'PosY', 'TargetClass', 'DistractorClass']
        mrk_pkt = nn.FillInstanceDataTable(fill_direction='bfill', grouping_fields=['TrialIdx'],
                                           fill_fields=trial_fill_fields)(data=mrk_pkt)

        # Load and parse the non-neural events.
        dat_filename = os.path.join(data_root, 'Cerebus_DATA',
                                    sess_info['name'], sess_info['date'], sess_info['nsx'])
        nsx_pkt = nn.ImportNSX(filename=dat_filename, load_signals=False,
                               load_spiketrains=False, load_waveforms=False)()
        nsx_pkt = misc.ParseDigitalWords()(data=nsx_pkt)
        mrk_pkt = nn.MergeStreams(data1=mrk_pkt, data2=nsx_pkt)()

        # For each non-neural event, find the corresponding marker, and save the NEV timestamp along with that marker.
        mrk_pkt = misc.GetAlignedNEVTimestamps()(data=mrk_pkt)

        # Recalculate trial outcome. Can use different gate-keeper events: 'Cue', 'Delay', 'Go'
        full_pkt = nn.MergeStreams(replace_if_exists=True)(data1=mrk_pkt, data2=eye_pkt, data3=gaze_events)
        full_pkt = misc.RecalculateOutcome(critical_event='Delay', performance_window=30)(data=full_pkt)

        # Save output
        nn.ExportH5(filename=os.path.join(preproc_root, sess_info['exp_code'] + '_behav.h5'))(data=full_pkt)
