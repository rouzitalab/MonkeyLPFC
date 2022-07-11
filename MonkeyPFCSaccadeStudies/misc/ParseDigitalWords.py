import logging
from datetime import datetime
import numpy as np
from neuropype.engine import *


logger = logging.getLogger(__name__)


class ParseDigitalWords(Node):
    # --- Input/output ports ---
    data = Port(None, Packet, "Data to process.", required=True,
                editable=False, mutating=True)

    field_name = StringPort(default="TargetValue")
    word_map = DictPort(default={
        'start_trial': 1,
        'stop_trial': 2,
        'lever_down': 3,
        'lever_release': 4,
        'flip_screen': 5,
        'start_my_stim': 6,
        'fixating': 7,
        'break_fixation': 8,
        'fix_req_on': 248,
        'fix_req_off': 249,
        'start_sentence': 250,
        'stop_sentence': 251,
        'start_recording': 252,
        'stop_recording': 253,
        'pause_recording': 254,
        'resume_recordings': 255
    })
    sentence_map = DictPort(default={'trial_ID': 1, 'file_ID': 2})

    @classmethod
    def description(cls):
        return Description(name='Parse Digital Words',
                           description="""
           Parse digital event values to get digital words and sentences.
                           """,
                           version='0.1',
                           license=Licenses.MIT)

    @data.setter
    def data(self, pkt):
        for n, c in enumerate_chunks(pkt, with_axes=(instance,), with_flags=Flags.is_event_stream):
            iax_dat = c.block.axes[instance].data
            words = iax_dat[self.field_name].astype(int)
            b_word_in_sentence = np.zeros_like(words, dtype=bool)

            # Identify sentences
            b_sent_starts = words == self.word_map['start_sentence']
            n_sentences = np.sum(b_sent_starts)

            # Output
            ev_times = np.nan * np.ones((n_sentences,))
            ev_markers = np.array([''] * n_sentences, dtype=object)
            ev_trial_id = np.nan * np.ones((n_sentences,))

            if np.any(b_sent_starts):
                sent_starts = np.where(b_sent_starts)[0]
                sent_stops = np.where(words == self.word_map['stop_sentence'])[0]
                # Cleanup unpaired starts/stops
                sent_starts = sent_starts[sent_starts < sent_stops[-1]]
                sent_stops = sent_stops[sent_stops > sent_starts[0]]
                # Prepare some of the easier variables
                sentence_type_lookup = {v: k for k, v in self.sentence_map.items()}
                ev_markers = np.array([sentence_type_lookup[words[_ + 1] ]for _ in sent_starts])
                ev_times = c.block.axes[instance].times[sent_starts]
                for s_ix, s_start in enumerate(sent_starts):
                    s_stop = sent_stops[s_ix]
                    b_word_in_sentence[s_start:s_stop+1] = True
                    sent_value = words[s_start+2:s_stop]
                    if ev_markers[s_ix] == 'file_ID':
                        supp = "".join([chr(_) for _ in sent_value])
                    elif ev_markers[s_ix] == 'trial_ID':
                        ev_trial_id[s_ix] = np.array([100, 1]) @ sent_value[-2:]  # trial_id
                        if sent_value[0] < 1980:
                            sent_value[0] += 2000
                        sent_value[5] = sent_value[5] % 60
                        dt = datetime(*sent_value[:6])
                        supp = str(dt)

            word_lookup = {v: k for k, v in self.word_map.items()}
            event_words = words[~b_word_in_sentence]
            ev_markers_ext = np.array([word_lookup[_] for _ in event_words])
            ev_times_ext = c.block.axes[instance].times[~b_word_in_sentence]
            ev_trial_id_ext = np.nan * np.ones_like(ev_times_ext)

            # Get trial_ID for start_trial event, which happens before corresponding trial_ID events.
            b_tr_start = ev_markers_ext == 'start_trial'
            new_trial_ids = ev_trial_id[np.searchsorted(ev_times, ev_times_ext[b_tr_start])]
            ev_trial_id_ext[b_tr_start] = new_trial_ids

            import pandas as pd
            df = pd.DataFrame({'Time': np.hstack((ev_times, ev_times_ext)),
                               'Marker': np.hstack((ev_markers, ev_markers_ext)),
                               'TrialID': np.hstack((ev_trial_id, ev_trial_id_ext))
                               })\
                .sort_values(by=['Time'])\
                .fillna(method='ffill', axis='index').dropna()\
                .infer_objects()
            df[['TrialID']] = df[['TrialID']].astype(np.int32)
            iax = InstanceAxis(df['Time'].values, data=df.drop(columns=['Time']).to_records(index=False))
            c.block = Block(data=np.nan * np.ones((len(iax),)), axes=(iax,))

        self._data = pkt
