import logging
import numpy as np
from neuropype.engine import *
from indl.metrics import dprime


logger = logging.getLogger(__name__)


class RecalculateOutcome(Node):
    # --- Input/output ports ---
    data = Port(None, Packet, "Data to process.", required=True,
                editable=False, mutating=True)

    critical_event = EnumPort(default="Go", domain=['Cue', 'Delay', 'Go'], help="""
        Trials with saccades that occur after the 'Target' event but before the critical_event are discarded.""")
    gaze_zero_recentered = BoolPort(default=True)
    performance_window = IntPort(default=30, help="Calculate performance metrics (dprime, bias, accuracy) "
                                                  "in a window of trials this size. First N-1 trials will have nan")

    @classmethod
    def description(cls):
        return Description(name='Recalculate Outcome',
                           description="""
                           Recalculate outcome of macaque pfc saccade trials
                           from (new) gaze events.
                           """,
                           version='0.1',
                           license=Licenses.MIT)

    @data.setter
    def data(self, pkt):
        ptb_n, ptb_chunk = find_first_chunk(pkt, with_flags=Flags.is_event_stream,
                                            name_startswith='markers')
        if ptb_n is not None:
            logger.info("Recalculating trial outcomes from behaviour data...")

            #  See Common_Matlab\Behaviour\getTrialBehavResult.m
            #  See Common_Matlab\Helper\getNewClass.m
            ptb_params = ptb_chunk.props['ptb_params']

            # Get the variables necessary to evaluate saccades based on gaze position.
            # gaze events units are degrees (rads?) of visual angle from the top-left of the screen.
            # Thus, to easily compare these units to our expectations, we need the
            # saccade start point (fixation), saccade end points (targets), and the
            # allowable radius for each in units of degrees of visual angle.
            # These variables, if given in ptbParams, are in units of pixels.

            # Conversion
            # Pixel size a.k.a. centimeters per pixel
            # pSize = .0702;  % This came from Florian's code.
            pSize = np.linalg.norm(ptb_params['subjectScreenSize'] / 10) / np.linalg.norm(ptb_params['subjectScreenResolution'])
            d2m = ptb_params['subjectScreenDistance'] / 10  # 100; # distance to monitor in centimeters
            # To convert any pixel measurement to degrees:
            # x_in_deg = atand(x_in_px * pSize / d2m);
            # x_in_px = d2m * tand(x_in_deg) / pSize;

            # Fixation point is in exact center of screen
            fixationPx = (ptb_params['subjectScreenResolution'] / 2).astype(int)  # From bottom/top? left, in px
            fixationDeg = np.rad2deg(np.arctan2(pSize * fixationPx, d2m))  # From bottom/top? left, in deg.v.a
            # fixationDeg = fixationPx/ptbParams.subjectScreenPixelPerDegree

            # Stimuli were squares with length 2*radius
            # fixation 'radius'
            fixationRadiusDeg = 6  # 4; %deg
            if 'fixPointfixRadius' in ptb_params:
                fixationRadiusDeg = np.rad2deg(np.arctan2(ptb_params['fixPointfixRadius'] * pSize, d2m))
                if type(fixationRadiusDeg) is np.ndarray:
                    fixationRadiusDeg = fixationRadiusDeg[-1]

            # target 'radius'
            targetRadiusDeg = 8  # deg
            if 'targetPointfixRadius' in ptb_params:
                targetRadiusDeg = np.rad2deg(np.arctan2(ptb_params['targetPointfixRadius'] * pSize, d2m))
                if type(targetRadiusDeg) is np.ndarray:
                    targetRadiusDeg = targetRadiusDeg[-1]

            # Process the event table to determine behaviour (saccade) outcomes.
            # This is better than using the ptb online-calculated trial outcomes because
            # saccade events were recalculated with better data(from eyelink),
            # we can relax saccade-timing requirements,
            # and we can accommodate some slight shifts in the eyeTracker.

            ev_dat = ptb_chunk.block.axes[instance].data
            ev_times = ptb_chunk.block.axes[instance].times

            trial_idx = ev_dat['TrialIdx']
            uq_tr_idx, inv_tr_idx = np.unique(trial_idx, return_inverse=True)
            n_trials = len(uq_tr_idx)
            n_evs = len(trial_idx)

            b_targs = ev_dat['Marker'] == 'Target'
            b_crits = ev_dat['Marker'] == self.critical_event
            b_sac = np.in1d(ev_dat['Marker'], ['Saccade', 'Gaze'])
            b_delay = ev_dat['Marker'] == 'Delay'
            b_go = ev_dat['Marker'] == 'Go'
            b_orig_occ = np.in1d(ev_dat['OutcomeCode'], [0, 9])

            # Recalculate the behavioural outcome and store in new fields.
            from collections import OrderedDict
            sac_ev = OrderedDict()  #{}
            # newOutcomeCode: -1 by default, 0 for success, 9 for saccades to distractor
            sac_ev['newOutcomeCode'] = -1 * np.ones((n_evs,), dtype=int)
            # critTime: The time of the critical event, generally used to indicate if a saccade was premature or not.
            sac_ev['critTime'] = np.nan * np.ones((n_evs,))
            # sacStartTime:
            sac_ev['sacStartTime'] = np.nan * np.ones((n_evs,))
            # sacClass:
            sac_ev['sacClass'] = -1 * np.ones((n_evs,), dtype=int)
            # runAcc
            sac_ev['runAcc'] = np.nan * np.ones((n_evs,))
            # runDPrime
            sac_ev['runDPrime'] = np.nan * np.ones((n_evs,))
            # runBias
            sac_ev['runBias'] = np.nan * np.ones((n_evs,))

            # Props for newOutcomeCode, critTime, sacStartTime, sacClass, runAcc, runDPrime, runBias
            sac_ev_props = [ValueProperty.CATEGORY & ValueProperty.INTEGER, ValueProperty.UNKNOWN,  # newOCC and crit
                            ValueProperty.UNKNOWN, ValueProperty.CATEGORY & ValueProperty.INTEGER,  # sac
                            ValueProperty.UNKNOWN, ValueProperty.UNKNOWN, ValueProperty.UNKNOWN]  # run

            for _ix, tr_idx in enumerate(uq_tr_idx):
                b_trial = ev_dat['TrialIdx'] == tr_idx

                # Mark trials with coinciding Delay and Go events as bad (IsGood = False) and skip processing.
                if np.any(np.in1d(ev_times[b_trial & b_delay], ev_times[b_trial & b_go])):
                    ev_dat['IsGood'][b_trial] = False
                    continue

                if np.any(b_trial & b_crits):
                    tr_sac_times = ev_times[b_trial & b_sac]
                    tr_crit_time = ev_times[b_trial & b_crits][0]
                    sac_ev['critTime'][b_trial] = tr_crit_time  # Fill all events for this trial with the critTime
                    tr_row = ev_dat[b_trial & b_crits][0]

                    # Saccades between target-onset and critical event are 'bad'
                    b_bad_sac = np.logical_and(tr_sac_times >= ev_times[b_trial & b_targs][0],
                                               tr_sac_times < tr_crit_time)
                    # Saccades between critical event and last event in trial are 'good'
                    b_good_sac = np.logical_and(tr_sac_times >= tr_crit_time,
                                                tr_sac_times <= ev_times[b_trial][-1])
                    # If we don't have any disqualifying saccades, and we have good saccades, calculate new outcome.
                    if ~np.any(b_bad_sac) and np.any(b_good_sac):
                        sac_rows = ev_dat[b_trial & b_sac][b_good_sac]
                        # (EndTime, Marker, Duration, Amp, PosX, StartX, PosY, StartY, TargetValue, IsGood)
                        start_pos_deg = np.vstack((sac_rows['StartX'], sac_rows['StartY']))
                        end_pos_deg = np.vstack((sac_rows['PosX'], sac_rows['PosY']))
                        if self.gaze_zero_recentered:
                            start_pos_deg = start_pos_deg + fixationDeg[:, None]
                            end_pos_deg = end_pos_deg + fixationDeg[:, None]

                        # Check which saccades started in the fixation point and finished in the target or distractor.
                        targ_deg = [np.rad2deg(np.arctan2(tr_row[_] * pSize, d2m)) for _ in ['TargetX', 'TargetY']]
                        b_fix_ok = [pointInRectangle(_, fixationDeg, fixationRadiusDeg) for _ in start_pos_deg.T]
                        b_in_targ = [pointInRectangle(_, targ_deg, 2*targetRadiusDeg) for _ in end_pos_deg.T]
                        dist_deg = [np.rad2deg(np.arctan2(tr_row[_] * pSize, d2m)) for _ in
                                    ['DistractorX', 'DistractorY']]
                        b_in_dist = [pointInRectangle(_, dist_deg, 2 * targetRadiusDeg) for _ in end_pos_deg.T]
                        b_targ_ok = np.any(b_in_targ) and not np.any(b_in_dist)
                        b_dist_ok = not np.any(b_in_targ) and np.any(b_in_dist)
                        if np.any(b_fix_ok) and (b_targ_ok or b_dist_ok):
                            sac_ev['newOutcomeCode'][b_trial] = 0 if b_targ_ok else 9
                            sac_ev['sacStartTime'][b_trial] = tr_sac_times[b_good_sac][np.where(b_targ_ok or b_dist_ok)[0][0]]
                            sac_ev['sacClass'][b_trial] = tr_row['TargetClass' if b_targ_ok else 'DistractorClass']
                        # elif np.any(np.in1d(sac_rows['OutcomeCode'], [0, 9])):
                        #     print("Why is OutcomeCode good but not my analysis?")

                    # If this trial has an old or new outcome code of 0 or 9, calculate its running behavioural perf.
                    # Use events that meet the following criteria: in this block, up to and including this trial,
                    # with crit event (1 per trial), original or new outcomecode in [0, 9].
                    b_eval_window = (ev_dat['Block'] == tr_row['Block']) & (ev_dat['TrialIdx'] <= tr_idx)\
                                    & b_crits & (b_orig_occ | np.in1d(sac_ev['newOutcomeCode'], [0, 9]))
                    if np.sum(b_eval_window) >= self.performance_window:
                        event_inds = np.where(b_eval_window)[0][-self.performance_window:]
                        targ_class = ev_dat['TargetClass'][event_inds]
                        targ_outcomes = sac_ev['newOutcomeCode'][event_inds]
                        b_use_old = ~np.in1d(targ_outcomes, [0, 9])
                        targ_outcomes[b_use_old] = ev_dat['OutcomeCode'][event_inds][b_use_old]
                        b_correct = targ_outcomes == 0  # targ_outcomes 0/9 --> 1/0
                        uq_targs, y_true = np.unique(targ_class, return_inverse=True)  #  0:3 / 4:7 --> 0 / 1
                        y_behav = y_true * b_correct + (1 - y_true) * ~b_correct
                        _dprime, bias, acc = dprime(y_true, y_behav)
                        sac_ev['runDPrime'][b_trial] = _dprime
                        sac_ev['runBias'][b_trial] = bias
                        sac_ev['runAcc'][b_trial] = acc

                    # if np.any(b_bad_sac):
                    #     print("{}: sacc during crit.".format(tr_idx))
                    # elif not b_fix_ok:
                    #     print("{}: bad fix".format(tr_idx))
                    # elif (not b_in_targ) and (not b_in_dist):
                    #     print("{}: bad sacc dest.".format(tr_idx))

            ptb_chunk.block.axes[instance].append_fields(list(sac_ev.keys()), list(sac_ev.values()),
                                                         props=sac_ev_props, overwrite_if_exists=True)
            self._data = pkt


def pointInRectangle(pt, rectCent, rectEdgeL):
    b_result = (rectCent[0] - rectEdgeL/2) <= pt[0] <= (rectCent[0] + rectEdgeL/2)
    b_result &= (rectCent[1] - rectEdgeL / 2) <= pt[1] <= (rectCent[1] + rectEdgeL / 2)
    return b_result
