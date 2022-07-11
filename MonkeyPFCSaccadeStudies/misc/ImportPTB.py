import logging
import scipy.io
import numpy as np
from neuropype.engine import *
from neuropype.utilities import cache
from neuropype.utilities.cloud import storage


logger = logging.getLogger(__name__)


class ImportPTB(Node):
    # --- Input/output ports ---
    filename = StringPort("", """Name of the event dataset.
                    """, is_filename=True, direction=IN)
    data = Port(None, Packet, "Output with markers chunk", required=True, direction=OUT)

    # options for cloud-hosted files
    cloud_host = EnumPort("Default", ["Default", "Azure", "S3", "Google",
                                      "Local", "None"], """Cloud storage host to
            use (if any). You can override this option to select from what kind of
            cloud storage service data should be downloaded. On some environments
            (e.g., on NeuroScale), the value Default will be map to the default
            storage provider on that environment.""")
    cloud_account = StringPort("", """Cloud account name on storage provider
            (use default if omitted). You can override this to choose a non-default
            account name for some storage provider (e.g., Azure or S3.). On some
            environments (e.g., on NeuroScale), this value will be
            default-initialized to your account.""")
    cloud_bucket = StringPort("", """Cloud bucket to read from (use default if
            omitted). This is the bucket or container on the cloud storage provider
            that the file would be read from. On some environments (e.g., on
            NeuroScale), this value will be default-initialized to a bucket
            that has been created for you.""")
    cloud_credentials = StringPort("", """Secure credential to access cloud data
            (use default if omitted). These are the security credentials (e.g.,
            password or access token) for the the cloud storage provider. On some
            environments (e.g., on NeuroScale), this value will be
            default-initialized to the right credentials for you.""")
    use_caching = BoolPort(False, """Enable caching.""", expert=True)

    @classmethod
    def description(cls):
        return Description(name='Import PyschoPhysicsToolbox data file'
                                ' from macaque saccade experiments',
                           description="""
                           """,
                           version='0.1',
                           license=Licenses.MIT)

    @Node.update.setter
    def update(self, v):
        record = cache.try_lookup(context=self, enabled=self.use_caching,
                                  verbose=True, data=self._data, state=None)
        if record.success():
            self._data = record.data
            return

        filename = storage.cloud_get(self.filename, host=self.cloud_host,
                                     account=self.cloud_account,
                                     bucket=self.cloud_bucket,
                                     credentials=self.cloud_credentials)

        logger.info("Loading Psychophysicstoolbox events from %s..." % filename)
        mat = scipy.io.loadmat(filename, squeeze_me=True)
        trial_data = mat['eData']['trial'][()]  # recarray

        # Do I need mat['DIO'].dtype.names? So far, no.

        # Normalize field names that vary across data sets.
        import numpy.lib.recfunctions as rfn
        trial_data = rfn.rename_fields(trial_data, {
            'trialStartTime': 'startTime',
            'trialStopTime': 'stopTime',
            'eyeSyncStartTime': 'eyeSyncTime',
            'trialID': 'ID',
            'timeofrelease': 'leverRelease',
            'cueTime': 'cuePresentedTime'
        })
        n_complete = np.sum(np.in1d(trial_data['outcomeCode'], [0, 9]))
        logger.info(f"Monkey completed {n_complete}/{len(trial_data)} trials ({100*n_complete / len(trial_data):.2f}%)")

        # Fix lack of stopTime on last trial --> use last flipScreen.
        if not isinstance(trial_data[-1]['stopTime'], float) or np.isnan(trial_data[-1]['stopTime']):
            trial_data[-1]['stopTime'] = max(trial_data[-1]['flipScreen'])

        # If needed, replace stopTime with trial duration
        if (trial_data[0]['stopTime'] - trial_data[0]['startTime']) > 0:
            trial_data['stopTime'] = trial_data['stopTime'] - trial_data['startTime']

        # I'm not sure about this next one...
        # In files with both M and A trial types,
        # M trials have targetChoice 1-8 (maybe 1-16) and class is empty.
        # A trials have targetChoice 0 or 1 and class is targetChoice+1.
        if len(np.unique(trial_data['expType'])) > 1:
            trial_data = rfn.drop_fields(trial_data, ['class'])
            trial_data = rfn.rename_fields(trial_data, {'targetChoice': 'class'})

        # Create a pandas dataframe with one row per each important event.
        # Each row within a trial will share many common values (see tr_map and check_map).
        # Rows within a trial will vary in the 'Marker', 'Time', and 'EyelinkTime' columns.
        import pandas as pd
        tr_map = {'ExperimentType': 'expType', 'Block': 'block', 'OutcomeCode': 'outcomeCode',
                  'TargetChoice': 'targetChoice', 'TrialIdx': 'expTrial', 'Class': 'class'}
        check_map = {'SR3errorStrategy': 'SR3errorStrategy', 'SR3trainingLevel': 'SR3trainingLevel',
                     'SR3InitialtrainingLevel': 'SR3InitialtrainingLevel'}
        ev_map = {'Start': 'startTime', 'Cue': 'cuePresentedTime', 'Delay': 'cueExtinguishedTime',
                  'Go': 'fixPointExtinguishedTime', 'Gaze': 'gazeEnteredTargetTime',
                  'Reward': 'rewardTime', 'Stop': 'stopTime'}

        col_names = ['Marker', 'Time', 'EyelinkTime'] + [k for k in tr_map.keys()] + [k for k in check_map.keys()]
        df = pd.DataFrame(columns=col_names)

        if 'SR3trainingLevel' in trial_data.dtype.names:
            tl_len = np.array([len(_) if not isinstance(_, (int, float)) else 1 for _ in trial_data['SR3trainingLevel']])
            if np.any(tl_len > 1):
                print("FIXME")  # TODO: handle trainingLevel as vector

        for tr_ix, tr_dat in enumerate(trial_data):
            # Collect trial-level variables
            tr = {k: tr_dat[v] for k, v in tr_map.items()}
            tr = {**tr, **{k: tr_dat[v] if (v in tr_dat.dtype.names and isinstance(tr_dat[v], (int, float)))
                                            else np.nan for k, v in check_map.items()}}

            trial_t0 = tr_dat['startTime']
            for ev_str, ev_v in ev_map.items():
                # Check why first trial Time is double but rest are OK.
                if isinstance(tr_dat[ev_v], float):  # Is there a faster way to check that it is not empty np.array?
                    ev_t_delta = tr_dat[ev_v] if ev_str != 'Start' else 0
                    ev_dict = {
                        'Marker': ev_str,
                        'Time': trial_t0 + ev_t_delta,
                        'EyelinkTime': tr_dat['eyeSyncTime']/1000. + ev_t_delta
                    }
                    df = df.append({**tr, **ev_dict}, ignore_index=True)

            # Also add two more events from flipScreen field that are not found in other trial fields.
            for fs_ix, fs_str in enumerate(['Fixation', 'Target']):
                if len(tr_dat['flipScreen']) > fs_ix:
                    ev_dict = {
                        'Marker': fs_str,
                        'Time': trial_t0 + tr_dat['flipScreen'][fs_ix],
                        'EyelinkTime': tr_dat['eyeSyncTime']/1000. + tr_dat['flipScreen'][fs_ix]
                    }
                    df = df.append({**tr, **ev_dict}, ignore_index=True)

        # Fixup dataframe: Sort by Time, attempt to auto-convert floats to ints.
        df = df.sort_values('Time').reset_index(drop=True).infer_objects()

        centre_point = mat['params']['subjectScreenResolution'][()] / 2

        # We also want, for each trial, to get the stimulus information.
        # This is a bit tricky because the relationship between PTB variables and the actual stimulus
        # depends on the experiment type. Some files have more than one experiment type.

        # Get per-trial indices into target_theta (targ/dist), cue_colour_map, and per trial radii
        # Cue Info. The PTB code uses trial classes for which the value indexes into a predefined colour-map, below.
        cue_colour_map = np.array(['r', 'g', 'b'])
        cue_matrix = np.array([[0, 1], [0, 2], [1, 0], [1, 2], [2, 0], [2, 1]])  # 1=red; 2=green; 3=blue.

        # Target Info for SR3
        # In Adam's code, the 1:8 part of the class is mapped to theta with
        # target_theta = 45*((1:8)-8-1); then the targetPt is calculated as
        # target_yx = centre_yx + radius*[-cosd(target_theta) sind(target_theta)];
        # Notice yx, not xy. Then he fliplr's the result.
        # This can be more simply represented the following way.
        target_theta = np.deg2rad([270, 315, 0, 45, 90, 135, 180, 225])
        target_str = np.array(['UU', 'UR', 'RR', 'DR', 'DD', 'DL', 'LL', 'UL'])  # The y-zero was up, so flip u/d.

        # Only needed for SR4: groups for general rules.
        context_options = np.array(['AnyUp', 'AnyDown', 'AnyLeft', 'AnyRight'])
        # [d;d;d],[u;u;u]
        # [l;l;l];[r;r;r]
        target_groups = np.array([[[7, 0, 1], [5, 4, 3]], [[7, 6, 5], [1, 2, 3]]])
        targdist_groups = np.stack((target_groups, np.flip(target_groups, axis=1)), axis=-1)

        # Additional columns added to the data frame:
        # NewType = 'M', 'A', 'SR3', or 'SR4'
        # CueColour = 'r', 'g', or 'b'
        # TrainingLevel = 0: full exp;              1: fix only;      2: fix & dist;
        #                 3: full exp & invis dist; 4: 3 + extra cue; 5: 4 + extra cue
        # TargetRule = 'UU', 'UR', 'RR', etc. OR 'AnyUp', 'AnyRight', etc.
        # Radius = target/distractor radius in polar coordinates.
        # TargetClass = in 1:8 for 8 locations, or in 1:16 with two annuli.
        # TargetPol = theta
        # Target[X,Y] = [x y] coordinates in screen pixels.
        # (also distClass, distPol, distXY)
        stim_info = {
            'NewType': np.array([''] * len(df), dtype=object),
            'CueColour': np.array(['_'] * len(df)),
            'TrainingLevel': np.nan * np.ones((len(df),)),
            'TargetRule': np.array([''] * len(df), dtype=object),
            'Radius': np.zeros((len(df)), dtype=int),
        }
        for td_str in ['Target', 'Distractor']:
            stim_info = {
                **stim_info,
                td_str + 'Class': np.zeros((len(df)), dtype=int),
                td_str + 'Theta': np.nan * np.ones((len(df),)),
                td_str + 'X': np.nan * np.ones((len(df),)),
                td_str + 'Y': np.nan * np.ones((len(df),)),
                td_str + 'Str': np.array([''] * len(df), dtype=object),
            }

        uq_exp, uq_idx = np.unique(df['ExperimentType'], return_inverse=True)
        for _ix, exp_type in enumerate(uq_exp):
            b_ev = uq_idx == _ix  # Trials in this experiment type.
            n_evs = np.sum(b_ev)

            # For most experiments, the trial type is constant throughout.
            trial_types = [{'M': 'M', 'A': 'A', 'C': 'SR3', 'D': 'SR4'}[exp_type]] * n_evs
            # But this isn't necessarily true for exp_type D.
            if exp_type == 'D':
                # Ultimately it comes down to whether the PTB code calls
                # dSR3T2c_setscreens2 or dSR4T2c_setscreens2. The code path is as
                # follows:
                # (start in saccadeGUImain)
                # if on probation -> doSR3Trial2e
                # if not on probation -> doSR3Trial2d
                #
                # Whether or not we are on probation can be changed throughout the
                # experiment with a key press -> saved to trials.SR3errorStrategy,
                # but only on the trial for which the key was pressed (how can we
                # know what the error strategy was at the beginning of the
                # experiment?!)
                #
                # doSR3Trial2e runs the same trial over and over until 2/3 are
                # successful, then it calls dSR3T2c_setscreens2
                # doSR3Trial2d runs dSR3T2c_setscreens2 if exp_type C, or
                # dSR4T2c_setscreen2 if exp_type D.
                #
                # so, exp_type C is always dSR3T2c_setscreens2
                if mat['params']['currSession'][()] <= 49:
                    # JerryLee sra3_2 before 090626 are actually SR3, but every sra3_2
                    # after that is SR4. I haven't yet spotted a sure-fire way to know
                    # that other than from the date.
                    trial_types = ['SR3'] * n_evs
                else:
                    # TODO: This block hasn't been checked with SR4 data.
                    print("FIXME: This block has not been checked.")
                    b_has_strat = ~df['SR3errorStrategy'].isna()
                    err_chng_to = df[b_has_strat]['SR3errorStrategy']
                    # err_chng_to = {err_chng_to.error}  # TODO: What's this?
                    # if we don't know the beginning error strategy, assume resample.
                    if ~b_has_strat[0]:
                        b_has_strat[0] = True
                        err_chng_to = ['resample'] + err_chng_to

                    is_prob = np.zeros((n_evs,), dtype=bool)
                    err_chng_id = np.where(b_has_strat)[0]
                    for e_ix in range(len(err_chng_id)):
                        if e_ix == len(err_chng_id):
                            this_slice = np.s_[err_chng_id[e_ix]:n_evs]
                        else:
                            this_slice = np.s_[err_chng_id[e_ix]:err_chng_id(e_ix + 1) - 1]

                        is_prob[this_slice] = err_chng_to[e_ix] == 'probation'

                    trial_types[is_prob] = ['SR3'] * sum(is_prob)

            if exp_type == 'M':
                radius = mat['params']['FMgridLength'] / 2*mat['params']['FMnumAnnuli:-1:1))']
                fix_jitter_xy = [0, 0]  # Add to fix/targ/distr XY
            elif exp_type == 'A':
                radius = mat['params']['targetDistance']
            elif exp_type == 'C':
                radius = [mat['params']['SR3radius'][()]]
                fix_jitter_xy = [mat['params'][_][()] for _ in ['SR3Xtranslate', 'SR3Ytranslate']]
            elif exp_type == 'D':
                radius = [mat['params']['SR3radius'][()][1]]
                fix_jitter_xy = [mat['params'][_][()][1] for _ in ['SR3Xtranslate', 'SR3Ytranslate']]

            # Training levels for each trial.
            # Determines whether target,distractor were used.
            training_levels = np.zeros((n_evs,))
            # TODO: Should SR3trainingLevel or SR3InitialtrainingLevel get priority here?
            b_tl_nan = df['SR3trainingLevel'].isna()
            training_levels[~b_tl_nan] = df['SR3trainingLevel'][~b_tl_nan]
            b_itl = b_tl_nan & ~df['SR3InitialtrainingLevel'].isna()
            training_levels[b_itl] = df['SR3InitialtrainingLevel'][b_itl]

            b_targ = training_levels != 1  # trial has target
            b_dist = ~np.in1d(training_levels, [1, 3, 4])  # trial has distractor
            # df[b_ev & ~b_dist].newType = 'CentreOut'

            #  method depends on trialType
            n_exp_trials = np.sum(b_ev)
            targ_ix = np.zeros((n_exp_trials,), dtype=int)
            dist_ix = np.zeros((n_exp_trials,), dtype=int)
            cue_colour_ix = np.zeros((n_exp_trials,), dtype=int)
            annulus_ix = np.zeros((n_exp_trials,), dtype=int)
            contexts = np.array([''] * n_exp_trials, dtype=object)

            classes = df[b_ev]['Class']

            uq_types, uq_inds = np.unique(trial_types, return_inverse=True)
            for type_ix, trial_type in enumerate(uq_types):
                b_type = uq_inds == type_ix
                type_classes = classes[b_type] - 1  # -1 because we use 0-based indexing here unlike Matlab 1-based.
                if trial_type == 'M':
                    flip_names = ['fixationOnset', 'targetOnset', 'fixationOffset', 'imperativeCue', 'saccadeEnd',
                                  'targetAcqLost']
                    # It was tough to follow the PTB code, I don't know if the above is correct.
                    # For "M" trials, we are looking for correct saccades typically after flipScreen 2 or 3
                    # ??? Class 0 is up, then around the face of a clock in 45
                    # degree increments until 7. If 8:15 present, those are the same except larger amplitude.
                    targ_ix[b_type] = type_classes % 8
                    annulus_ix[b_type] = np.floor(type_classes / 8).astype(int)
                    cue_colour_ix[b_type] = 0
                elif trial_type == 'A':
                    logger.warning("Not yet implemented exp_type A.")
                elif trial_type == 'SR3':
                    flip_names = ['fixationOnset', 'targetOnset', 'cueOnset', 'cueOffset', 'fixationOffset',
                                  'targetAcqLost']  # fixationOffset is the imperative cue.
                    # trial class = target_config + 4*cue_index + 8*cue_row
                    # target_config is in 0:3, representing 4 direction pairs, rand per block
                    # I rename target_config -> theta_targ_ix
                    # cue_index is in 0:1, indexing which half target is in.
                    # --> targets are in [270 315 0 45 90 135 180 225]
                    # cue_row is in 0:length(cue_matrix)-1, chosen rand at the start of each block
                    # targetChoice is in 0:7, = 4*cue_index + target_config
                    cue_row = np.floor(type_classes / 8).astype(int)
                    cue_index = np.floor((type_classes - 8 * cue_row) / 4).astype(int)
                    target_config = type_classes - 8 * cue_row - 4 * cue_index
                    targ_ix[b_type] = 4 * cue_index + target_config
                    dist_ix[b_type] = targ_ix[b_type] + 4 * (targ_ix[b_type] < 4) - 4 * (targ_ix[b_type] > 3)
                    cue_colour_ix[b_type] = np.ravel_multi_index((cue_row, cue_index), cue_matrix.shape)
                    annulus_ix[b_type] = 0
                    contexts[b_type] = target_str[targ_ix[b_type].tolist()]
                elif trial_type == 'SR4':
                    flip_names = ['fixationOnset', 'targetOnset', 'cueOnset', 'cueOffset', 'fixationOffset',
                                  'saccadeEnd', 'targetAcqLost']

                    # trial class = 108*(r-1) + 54*(y-1) + 9*(x-1) + z

                    # r is cueColumn; 1 for u|l correct, 2 for d|r corr, rand per trial
                    r = np.ceil(type_classes / 108)
                    type_classes = type_classes - 108 * (r - 1)

                    # y is 1 for u/d, 2 for l/r, chosen randomly per block
                    y = np.ceil(type_classes/54)
                    temp = type_classes - 54*(y-1)

                    # x is cue_row in cue_matrix, in 1:6, (determines targ/dist colour pairings per block)
                    x = np.ceil(type_classes / 9)
                    type_classes = type_classes - 9 * (x - 1)

                    # z in 1:9, rand per block. Indexes into target,distractor pairings
                    # t order: 1 2 3 1 2 3 1 2 3 into targdist_groups(t_ix, r, y, 1)
                    # d order: 1 1 1 2 2 2 3 3 3 into targdist_groups(d_ix, r, y, 2)
                    z = type_classes
                    t_ix = z % 3
                    t_ix[t_ix == 0] = 3
                    d_ix = np.ceil(z/3)

                    targ_ix[b_type] = targdist_groups[
                        np.ravel_multi_index((t_ix, r, y, np.ones(t_ix.shape)), targdist_groups.shape)]
                    dist_ix[b_type] = targdist_groups[
                        np.ravel_multi_index((d_ix, r, y, 2 * np.ones(t_ix.shape)), targdist_groups.shape)]
                    cue_colour_ix[b_type] = np.ravel_multi_index((x, r), cue_matrix.shape)
                    annulus_ix[b_type] = 0
                    contexts[b_type] = context_options[(y-1)*2 + r]
                else:
                    logger.error("trial type not recognized.")

            # Start filling in stim_info
            stim_info['NewType'][b_ev] = trial_types
            stim_info['CueColour'][b_ev] = cue_colour_map[cue_matrix.flatten()[cue_colour_ix]]
            stim_info['TrainingLevel'][b_ev] = training_levels
            stim_info['TargetRule'][b_ev] = contexts
            stim_info['Radius'][b_ev] = np.array(radius)[annulus_ix]

            # Process fields for targets and distractors.
            targdist_ix = np.stack((targ_ix, dist_ix), axis=-1)
            b_td = np.stack((b_targ, b_dist), axis=-1)
            for td, tdstr in enumerate(['Target', 'Distractor']):
                b_is_td = b_td[:, td]         # If trial had a target/distractor
                this_ix = targdist_ix[b_is_td, td]  # Indices into target_theta

                # Stimulus class of 8 (or 16) possible locations.
                tr_class = np.zeros((n_evs,), dtype=int)
                tr_class[b_is_td] = this_ix
                tr_class[annulus_ix > 0] = tr_class[annulus_ix > 0] * (1 + annulus_ix[annulus_ix > 0])
                stim_info[tdstr + 'Class'][b_ev] = tr_class

                # Polar coordinates: theta and radius. Radius does not depend on t/d and is done already.
                tr_theta = np.nan * np.ones((n_evs,))
                tr_theta[b_is_td] = target_theta[this_ix]
                stim_info[tdstr + 'Theta'][b_ev] = tr_theta

                # Coordinates in screen pixels
                stim_info[tdstr + 'X'][b_ev] = stim_info['Radius'] * np.cos(tr_theta)\
                                                  + centre_point[0] + fix_jitter_xy[0]
                stim_info[tdstr + 'Y'][b_ev] = stim_info['Radius'] * np.sin(tr_theta)\
                                                  + centre_point[1] + fix_jitter_xy[1]

                # String for direction
                tr_str = np.array([''] * n_evs, dtype=object)
                tr_str[b_is_td] = target_str[this_ix]
                stim_info[tdstr + 'Str'][b_ev] = tr_str

        df = df.assign(**stim_info).infer_objects()

        # Save online eyetracker calibration adjustments
        if 'CALADJ' in mat:
            gaze_calib_adjust = mat['CALADJ']
        elif 'eyePosCalibAdjust' in mat['eData'].dtype.names:
            gaze_calib_adjust = mat['eData']['eyePosCalibAdjust'][()]
        else:
            gaze_calib_adjust = None

        # Save the output in Neuropype format.
        iax = InstanceAxis(df['Time'].values, data=df.drop(columns=['Time']).to_records(index=False))
        ev_blk = Block(data=np.nan * np.ones((len(iax),)), axes=(iax,))
        ev_props = {Flags.is_event_stream: True,
                    'gaze_calib_adjust': (gaze_calib_adjust,),  # Need to hide in tuple so it doesn't break enumerator
                    'ptb_params': {**{k: mat['params'][k][()] for k in mat['params'].dtype.names},
                                   'flip_names': flip_names}
                    }
        self._data = Packet(chunks={'markers': Chunk(block=ev_blk, props=ev_props)})

        record.writeback(data=self._data)
