import numpy as np
__all__ = ['get_data_root', 'get_chan_list', 'load_macaque_pfc', 'sess_infos']


TARG_ANGLES = np.asarray([2, 1, 0, 7, 6, 5, 4, 3])  # * np.pi/4
TARG_MAP = {'Target-' + str(_ + 1): TARG_ANGLES[_] for _ in range(8)}
sess_infos = [
    {'name': 'JerryLee', 'bank': 'A', 'name_short': 'j', 'date': '090601', 'exp_code': 'sra3_2_j_037_00+03', 'nsx': 'datafile003.ns2'},
    # {'name': 'JerryLee', 'bank': 'C', 'name_short': 'j', 'date': '090622', 'exp_code': 'sra3_2_j_049_00+02', 'nsx': 'datafile002.ns2'},
    {'name': 'JerryLee', 'bank': 'A', 'name_short': 'j', 'date': '090623', 'exp_code': 'sra3_1_j_050_00+', 'nsx': 'datafile002.ns2'},
    {'name': 'JerryLee', 'bank': 'B', 'name_short': 'j', 'date': '090624', 'exp_code': 'sra3_1_j_051_00+', 'nsx': 'datafile002.ns2'},
    {'name': 'JerryLee', 'bank': 'C', 'name_short': 'j', 'date': '090625', 'exp_code': 'sra3_1_j_052_00+', 'nsx': 'datafile003.ns2'},
    # {'name': 'Marty', 'bank': 'A', 'name_short': 'm', 'date': '090908', 'exp_code': 'sra3_1_m_070_00+02', 'nsx': 'datafile001.ns2'},
    # {'name': 'Marty', 'bank': 'A', 'name_short': 'm', 'date': '090914', 'exp_code': 'sra3_1_m_074_00+01', 'nsx': 'datafile001.ns2'},
    # {'name': 'Marty', 'bank': 'B', 'name_short': 'm', 'date': '090918', 'exp_code': 'sra3_1_m_076_00+01', 'nsx': 'datafile001.ns2'},
    {'name': 'Marty', 'bank': 'C', 'name_short': 'm', 'date': '090919', 'exp_code': 'sra3_1_m_077_00+01', 'nsx': 'datafile001.ns2'},
    # {'name': 'Marty', 'bank': 'B', 'name_short': 'm', 'date': '090920', 'exp_code': 'sra3_1_m_078_00+01', 'nsx': 'datafile001.ns2'},  # Not enough correct trials
    # {'name': 'Marty', 'bank': 'C', 'name_short': 'm', 'date': '090921', 'exp_code': 'sra3_1_m_079_00+01', 'nsx': 'datafile001.ns2'},
    # {'name': 'Marty', 'bank': 'A', 'name_short': 'm', 'date': '090922', 'exp_code': 'sra3_1_m_080_00+01', 'nsx': 'datafile001.ns2'},  # Not enough correct trials
    {'name': 'Marty', 'bank': 'B', 'name_short': 'm', 'date': '090925', 'exp_code': 'sra3_1_m_081_00+01', 'nsx': 'datafile001.ns2'},
    {'name': 'Marty', 'bank': 'A', 'name_short': 'm', 'date': '090926', 'exp_code': 'sra3_1_m_082_00+01', 'nsx': 'datafile001.ns2'},
    {'name': 'Marty', 'bank': 'C', 'name_short': 'm', 'date': '090927', 'exp_code': 'sra3_1_m_083_00+01', 'nsx': 'datafile001.ns2'},
    # {'name': 'Marty', 'bank': 'A', 'name_short': 'm', 'date': '090928', 'exp_code': 'sra3_1_m_084_00+01', 'nsx': 'datafile001.ns2'},
    # {'name': 'Marty', 'bank': 'B', 'name_short': 'm', 'date': '090929', 'exp_code': 'sra3_1_m_085_00+01', 'nsx': 'datafile001.ns2'}
]
# D:\SachsLab\Studies\Monkey\MonkeyStudies\SessionInfo\ArrayChannelMaps
map_files = {
    'JerryLee': 'JerryLee-1025-0359.cmp',
    'Marty': 'Marty-1025-0362.cmp'
}


def get_data_root():
    import os

    # TODO: Add different data roots for different computers.
    check_dirs = [
        # SachsLab Tower - Windows
        os.path.abspath(r'E:\Chad\The Ottawa Hospital\Sachs Lab - Documents\RawData\Monkey'),
        # SachsLab Tower - Linux
        os.path.abspath('/media/SachsLab/Chad/The Ottawa Hospital/Sachs Lab - Documents/RawData/Monkey'),
        # Chad:
        os.path.abspath(r'G:\SachsLab\RawData\Common-Monkey'),  # Chad home - Windows
        os.path.abspath('/media/chad/STORE/SachsLab/RawData/Common-Monkey'),  # Chad home - Linux
        os.path.abspath('/Volumes/CHADPORTABLE/SachsLab/Data/Common-Monkey')  # Chad Macbook w/ external HDD
    ]

    for datadir in check_dirs:
        if os.path.isdir(datadir):
            output = datadir
            break
    else:
        from tkinter import filedialog
        output = filedialog.askopenfilename(initialdir=os.path.dirname(os.getcwd()))
    return output


def get_chan_list(n_chans=32):
    import os
    import numpy as np
    if n_chans == 32:
        fname = 'chan_list_32.txt'
    elif n_chans == 64:
        fname = 'chan_list_64.txt'
    else:
        raise ValueError("chan list unavailable for provided n_chans")
    fname = os.path.join(os.path.dirname(__file__), fname)
    return np.loadtxt(fname, dtype='U')


def load_macaque_pfc(data_path, sess_id, x_chunk='spikerates', zscore=False,
                     valid_outcomes=(0,), dprime_range=(-np.inf, np.inf),
                     time_range=(-np.inf, np.inf),
                     min_block_length=0,
                     y_type=False,
                     sparse_to_bool=True,
                     resample_X=1,
                     samples_last=False,
                     verbose=0):
    """
    Load a file from the macaque pfc saccade dataset.
    :param data_path: Path object pointing to the root of the data directory
    :param sess_id: e.g., sra3_1_m_074_0001
    :param x_chunk: Type of data to return. Options are 'analogsignals' (i.e. LFPs), 'spikerates', 'spiketrains', 'gaze'
    :param zscore: Set True to center and standardize X data per-channel. Not available with spiketrains.
    :param valid_outcomes: tuple, should contain 0 (correct trials), and/or 9 (saccades to wrong target)
    :param dprime_range: 2-tuple, limits of acceptable dprime values. Set the first value to -np.inf to allow np.nan.
                         70% for both classes is dprime=1.0488.
    :param time_range: 2-tuple of time range to further limit data, in seconds.
    :param y_type: Data to return in y. Possible values:
        'pair and choice' for np.array of (target_pair, choice_within_pair)
        'encoded input' for ...
        any other string will be interpreted as a key into the instance axis data and will return that data.
        If evaluates to False (0, None, False), then default 'sacClass' is used.
    :param sparse_to_bool: If data is sparse, convert it to boolean events and discard the values associated with the ev
    :param resample_X: integer factor to resample X at. Uses np.nanmean
    :param samples_last: Set to True to transpose the time axis with the channels axis.
    :param verbose:
    :return: X, Y, ax_info
    """
    from indl.fileio import from_neuropype_h5
    sess_id = sess_id.replace("+", "")
    test_file = data_path / (sess_id + '_segmented.h5')
    chunks = from_neuropype_h5(test_file, chunk_names=[x_chunk])
    chunk_names = [_[0] for _ in chunks]
    chunk = chunks[chunk_names.index(x_chunk)][1]
    ax_types = [_['type'] for _ in chunk['axes']]
    inst_ix, time_ix, space_ix = (ax_types.index(_) for _ in ('instance', 'time', 'space'))
    instance_axis = chunk['axes'][inst_ix]

    # Keep trials based on dprime and outcome codes
    b_keep_inst = np.ones((len(instance_axis['data'].index),), dtype=bool)
    if 'runDPrime' in instance_axis['data'] and (dprime_range[0] > -np.inf or dprime_range[1] < np.inf):
        dprimes = instance_axis['data']['runDPrime'].values
        np.warnings.filterwarnings('ignore')
        b_dprime = np.logical_and(dprimes >= dprime_range[0], dprimes <= dprime_range[1])
        if dprime_range[0] == -np.inf:  # allow nan
            b_dprime = np.logical_or(b_dprime, np.isnan(dprimes))
        b_keep_inst = np.logical_and(b_keep_inst, b_dprime)

    if 'newOutcomeCode' in instance_axis['data']:
        occs = instance_axis['data']['newOutcomeCode'].values
        b_keep_inst = np.logical_and(b_keep_inst, np.in1d(occs, valid_outcomes))

    if min_block_length > 0:
        block_idx = instance_axis['data']['Block']
        block_switches = np.hstack(([True], np.diff(block_idx) > 0))
        block_switch_inds = np.where(block_switches)[0]
        trials_in_block = np.zeros((len(block_idx),))
        for bl_ix, sw_ind in enumerate(block_switch_inds):
            next_blk_start = block_switch_inds[bl_ix + 1] if bl_ix < (len(block_switch_inds) - 1) else len(block_idx)
            trials_in_block[sw_ind:] = next_blk_start - sw_ind
        b_keep_inst = np.logical_and(b_keep_inst, trials_in_block > min_block_length)

    tvec = chunk['axes'][time_ix]['times']
    b_tvec = np.logical_and(tvec >= time_range[0], tvec <= time_range[1])

    X = chunk['data'][b_keep_inst, :, :]
    X = np.nan_to_num(X)
    X = np.take(X, np.where(b_tvec)[0], axis=time_ix)

    ax_info = {'instance_data': instance_axis['data'][b_keep_inst],
               'instance_times': instance_axis['times'][b_keep_inst],
               'fs': chunk['axes'][time_ix]['nominal_rate'],
               'timestamps': chunk['axes'][time_ix]['times'][b_tvec],
               'channel_names': chunk['axes'][space_ix]['names'],
               'channel_locs': chunk['axes'][space_ix]['positions']
               }

    if time_ix > space_ix:
        # LFPs were stored with axes order reversed.
        X = np.transpose(X, (0, 2, 1))

    if sparse_to_bool and x_chunk == 'spiketrains':
        X = X != 0

    if zscore:
        if x_chunk == 'spiketrains':
            raise ValueError("Cannot z-score boolean spiketrains.")
        # Centre and standardize across all trials*samples
        tmp = np.reshape(X, (X.shape[0] * X.shape[1], X.shape[2]))
        X = (X - np.mean(tmp, axis=0)) / np.std(tmp, axis=0)

    if resample_X > 1:
        n_pad = resample_X - (X.shape[1] % resample_X)
        X = np.concatenate((X, np.nan*np.ones((X.shape[0], n_pad, X.shape[2]))), axis=1)
        X = np.reshape(X, (X.shape[0], X.shape[1] // resample_X, resample_X, X.shape[2]))
        X = np.nanmean(X, axis=-2)

        tvec = ax_info['timestamps']
        tvec = np.concatenate((tvec, np.nan*np.ones((n_pad,))))
        tvec = np.nanmean(np.reshape(tvec, (-1, resample_X)), axis=-1)
        ax_info['timestamps'] = tvec
        ax_info['fs'] = ax_info['fs'] / resample_X

    if y_type == 'pair and choice':
        Y = ax_info['instance_data']['sacClass'].values.reshape(-1, 1)
        Y = np.concatenate((Y % 4, Y > 3), axis=1).astype(int)
    elif y_type == 'encoded input':
        # 11-channel representation of the visual stimuli.
        n_trials = len(ax_info['instance_data'])
        tvec = ax_info['timestamps']
        n_times = len(tvec)
        Y = np.zeros((n_trials, n_times, 11))

        ch_offset = 0

        # First 4 channels give the 1-hot encoded target pair (sacClass % 4), but only when the target was presented
        # assuming target onset at t=0 until the end of the trial.
        b_t_targ = tvec >= 0
        for tr_ix, targ_pair in enumerate(ax_info['instance_data']['sacClass'].values % 4):
            Y[tr_ix, b_t_targ, ch_offset + targ_pair] = 1
        # assert np.all(np.any(np.any(Y[:, :, :4], axis=-1), axis=-1))
        ch_offset += 4

        # Next 3 channels give the 1-hot encoded colour cue (r, g, b)
        # assuming colour cue onset at t=# until t=#
        b_t_col = np.logical_and(tvec >= 0.25, tvec <= 1.25)
        cue_colours = ax_info['instance_data']['CueColour'].values
        for tr_ix, cc in enumerate(cue_colours):
            Y[tr_ix, b_t_col, ch_offset + ['r', 'g', 'b'].index(cc)] = 1
        # assert np.all(np.any(np.any(Y[:, :, 4:7], axis=-1), axis=-1))
        ch_offset += 3

        # Next 3 channels give the 1-hot encoded "rule" which is constant through the block.
        # This indicates which colour directed the animal to respond to the "preferred" target, where the preferred
        # targets were one of (TargetStr: TargetClass): UU:0, UR:1, RR:2, UL:7
        blk_ids = ax_info['instance_data']['Block'].values
        b_pref = np.in1d(ax_info['instance_data']['TargetStr'].values, ['UU', 'UR', 'RR', 'UL'])
        for blk_ix, blk_idx in enumerate(np.unique(blk_ids)):
            b_tr_blk = blk_ids == blk_idx
            b_blk_pref_colour = np.logical_and(b_pref, b_tr_blk)
            if np.any(b_blk_pref_colour):
                pref_cue_col = cue_colours[b_blk_pref_colour][0]
            else:
                # No trials in preferred targets so rule is unknown.
                # Try to find out from other trials in block that were removed.
                b_x_blk = instance_axis['data']['Block'].values == blk_idx
                b_x_pref = np.in1d(instance_axis['data']['TargetStr'].values, ['UU', 'UR', 'RR', 'UL'])
                b_blk_pref_colour = np.logical_and(b_x_pref, b_x_blk)
                if np.any(b_blk_pref_colour):
                    pref_cue_col = instance_axis['data']['CueColour'].values[b_blk_pref_colour][0]
                else:
                    print("TODO: Guess")
            Y[b_tr_blk, :, ch_offset + ['r', 'g', 'b'].index(pref_cue_col)] = 1
        # assert np.all(np.any(np.any(Y[:, :, 7:10], axis=-1), axis=-1))
        ch_offset += 3

        # Next, the fixation point. This acts as a 'hold' signal.
        # This was on from the beginning of the trial until 300 msec after colour offset, or +1.55
        b_t_fix = tvec < 1.55
        Y[:, b_t_fix, ch_offset] = 1
        ch_offset += 1

    else:
        if not y_type or not isinstance(y_type, str):
            y_type = 'sacClass'
        Y = ax_info['instance_data'][y_type].values.reshape(-1, 1)

    if verbose:
        print(f"Found {len(ax_info['instance_data'])} trials, "
              f"{len(ax_info['timestamps'])} timestamps"
              f"({ax_info['timestamps'][0]} to {ax_info['timestamps'][-1]} at {ax_info['fs']} Hz), "
              f"{X.shape[-1]} channels")
        print(f"Returning Y as {y_type} with shape {Y.shape}.")
        print(f"Axis info has: {ax_info.keys()}")

    if samples_last:
        X = np.transpose(X, (0, 2, 1))
        if y_type == 'encoded input':
            Y = np.transpose(Y, (0, 2, 1))

    return X, Y, ax_info


def dec_from_enc(encoded_stims):
    # We create a decision output for each of the 8 channels.
    # This is an output time series that reflects the available knowledge of what the output target is.
    # After target onset, the two possible targets have a value of 0.5, and after the color cue comes on
    # the correct target has a value of 1.0 and all others are 0.
    p_targ = np.zeros((encoded_stims.shape[0], encoded_stims.shape[1], 8))
    for tr_ix, tr_dat in enumerate(encoded_stims):
        targ_pair = tr_dat[:, :4]
        cue_col = tr_dat[:, 4:7]
        rule = tr_dat[:, 7:10]
        # fix = tr_dat[:, 10]
        # Create the output
        # First, from the targets.
        pair_ix = np.where(targ_pair)[1][0]
        b_targ_on = np.any(targ_pair, axis=1)
        p_targ[tr_ix, b_targ_on, pair_ix] = 0.5
        p_targ[tr_ix, b_targ_on, pair_ix + 4] = 0.5
        # Second, from the colour cue
        cue_onset = np.where(cue_col)[0][0]
        cue_ix = np.where(cue_col)[1][0]
        rule_ix = np.where(rule)[1][0]
        # If rule_ix and cue_ix match, then the target is one of UU:0, UR:1, RR:2, UL:7
        is_anti = (pair_ix in [3]) if rule_ix == cue_ix else (pair_ix in [0, 1, 2])
        targ_ix = pair_ix + 4 * is_anti
        p_targ[tr_ix, cue_onset:, :] = 0
        p_targ[tr_ix, cue_onset:, targ_ix] = 1
    return p_targ

if __name__ == '__main__':
    from pathlib import Path
    if Path.cwd().stem == 'misc':
        import os
        os.chdir('..')
    datadir = Path.cwd() / 'StudyLocationRule' / 'Data' / 'Preprocessed'
    res_tuple = load_macaque_pfc(datadir, 'sra3_1_j_050_00+',
                                 x_chunk='spikerates',
                                 valid_outcomes=(0,),
                                 zscore=True,
                                 dprime_range=(1.0, np.inf),
                                 time_range=(-np.inf, 1.45),
                                 verbose=True,
                                 min_block_length=10,
                                 y_type='encoded input',
                                 resample_X=5)
    print(res_tuple[2])  # X_ax_info
    knowledge = dec_from_enc(res_tuple[1])
    print(knowledge)
