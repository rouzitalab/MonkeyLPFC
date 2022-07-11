import logging
from datetime import datetime
import numpy as np
from neuropype.engine import *


logger = logging.getLogger(__name__)


class GetAlignedNEVTimestamps(Node):
    # --- Input/output ports ---
    data = Port(None, Packet, "Data to process.", required=True,
                editable=False, mutating=True)

    @classmethod
    def description(cls):
        return Description(name='Get Aligned NEV Timestamps',
                           description="""
                           """,
                           version='0.1',
                           license=Licenses.MIT)

    @data.setter
    def data(self, pkt):
        mrk_n, mrk_c = find_first_chunk(pkt, name_startswith='markers')
        nev_n, nev_c = find_first_chunk(pkt, name_startswith='events')

        if mrk_n is not None and nev_n is not None:
            nev_iax = nev_c.block.axes[instance]
            mrk_iax = mrk_c.block.axes[instance]

            # In some recordings a new PTB trial is started but the user interrupts it
            # and the partial trial information is not recorded in the nev though it is recorded in the ptb.
            while mrk_iax.data['TrialIdx'][-1] > nev_iax.data['TrialID'][-1]:
                last_trial_idx = mrk_iax.data['TrialIdx'][-1]
                b_keep = mrk_iax.data['TrialIdx'] < last_trial_idx
                mrk_c.block = mrk_c.block[..., instance[b_keep], ...]
                mrk_iax = mrk_c.block.axes[instance]

            # NEV timestamps will be stored in nev_timestamps, and eventually added as a column to marker inst axis.
            nev_timestamps = np.nan * np.ones_like(mrk_iax.times)

            # Align trial starts.
            b_nev_starts = nev_iax.data['Marker'] == 'start_trial'
            b_mrk_starts = mrk_iax.data['Marker'] == 'Start'
            nev_timestamps[b_mrk_starts] = nev_iax.times[b_nev_starts]

            # Align known flip screen events.
            nev_flips = mrk_c.props['ptb_params']['flip_names']
            flip_mrk_map = {
                'targetOnset': 'Target',
                'cueOnset': 'Cue',
                'cueOffset': 'Delay',
                'fixationOffset': 'Go'
            }
            b_nev_flip = nev_iax.data['Marker'] == 'flip_screen'
            for _ix, tr_idx in enumerate(np.unique(mrk_iax.data['TrialIdx'])):
                b_mrk_tr = mrk_iax.data['TrialIdx'] == tr_idx
                b_nev_flip_tr = b_nev_flip & (nev_iax.data['TrialID'] == tr_idx)
                for flip_ix in range(np.sum(b_nev_flip_tr)):
                    flip_name = nev_flips[flip_ix]
                    if flip_name in flip_mrk_map:
                        b_mrk_flip = b_mrk_tr & (mrk_iax.data['Marker'] == flip_mrk_map[flip_name])
                        nev_timestamps[b_mrk_flip] = nev_iax.times[b_nev_flip_tr][flip_ix]

            if False:
                import matplotlib.pyplot as plt
                plt.scatter(nev_timestamps[~np.isnan(nev_timestamps)], mrk_iax.times[~np.isnan(nev_timestamps)])
                plt.show()

            mrk_c.block.axes[instance].append_fields(('NEVTimestamps',), (nev_timestamps,), overwrite_if_exists=True)
            pkt.chunks[mrk_n] = mrk_c

        self._data = pkt
