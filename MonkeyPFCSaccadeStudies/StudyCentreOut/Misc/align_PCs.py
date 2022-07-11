from pathlib import Path
import scipy.io
import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt  # For debug plotting


fname = 'ForKaiM3PCVecsToBeAligned.mat'
ref_chan = None  # Set to None to infer reference channel. Otherwise `17` is good.
pc_cutoff = 10
f_vec = np.array([
    1.953125, 3.90625, 5.859375, 7.8125, 9.765625, 11.71875, 13.671875, 15.625, 17.578125, 19.53125,
    21.484375, 23.4375, 25.390625, 27.34375, 29.296875, 31.25, 33.203125, 35.15625, 37.109375, 39.0625,
    41.015625, 42.96875, 44.921875, 46.875, 48.828125, 50.78125, 52.734375, 54.6875, 56.640625, 64.453125,
    66.40625, 68.359375, 70.3125, 72.265625, 74.21875, 76.171875, 78.125, 80.078125, 82.03125, 83.984375,
    85.9375, 87.890625, 89.84375, 91.796875, 93.75, 95.703125, 97.65625, 99.609375, 101.5625, 103.515625,
    105.46875, 107.421875, 109.375, 111.328125, 113.28125, 115.234375, 123.046875, 125., 126.953125, 128.90625,
    130.859375, 132.8125, 134.765625, 136.71875, 138.671875, 140.625, 142.578125, 144.53125, 146.484375, 148.4375,
    150.390625, 152.34375, 154.296875, 156.25, 158.203125, 160.15625, 162.109375, 164.0625, 166.015625, 167.96875,
    169.921875, 171.875, 173.828125, 175.78125, 183.59375, 185.546875, 187.5, 189.453125, 191.40625, 193.359375,
    195.3125, 197.265625, 199.21875])


def debug_plot(pc_weights, n_plots=4, color_ignore=5):
    fig, axes = plt.subplots(n_plots, 1, figsize=(8, 11))
    for pc_ix, ax in enumerate(axes):
        for_clim = pc_weights[pc_ix, :, color_ignore:]
        ax.imshow(pc_weights[pc_ix, :, :], interpolation='none', vmin=np.min(for_clim), vmax=np.max(for_clim),
                  extent=[f_vec[0], f_vec[-1], pc_weights.shape[1], 1], aspect=2)
    plt.show()


if __name__ == "__main__":
    path = Path(__file__).parents[1] / 'Data' / 'Preprocessed' / fname
    mat_dict = scipy.io.loadmat(path, squeeze_me=True)

    pc_vecs = mat_dict['pc_vecs']
    n_pcs, n_chans, _ = pc_vecs.shape

    debug_plot(pc_vecs)

    if ref_chan is None:
        # Find which channel has minimum cost when comparing its abs weights to the average abs weights
        avg_abs = np.mean(np.abs(pc_vecs), axis=1)
        ref_cost = np.inf
        for ch_ix in range(n_chans):
            test_ref = np.abs(pc_vecs[:, ch_ix])
            sim = np.dot(avg_abs, test_ref.T)
            cost = 1 - np.abs(sim)
            cost = np.sum(np.diag(cost))
            if cost < ref_cost:
                ref_chan = ch_ix
                ref_cost = cost

    ref_mat = pc_vecs[:pc_cutoff, ref_chan, :]
    ref_mat = ref_mat / np.linalg.norm(ref_mat)

    align_out = np.nan * np.ones((pc_cutoff, n_chans, n_pcs))
    for chan_ix in range(n_chans):
        test_mat = pc_vecs[:pc_cutoff, chan_ix, :]
        test_mat = test_mat / np.linalg.norm(test_mat)
        sim = np.dot(ref_mat, test_mat.T)
        cost = 1 - np.abs(sim)
        row_inds, col_inds = scipy.optimize.linear_sum_assignment(cost)
        if not np.array_equal(row_inds, np.arange(len(row_inds)).astype(int)):
            print("oops")
        align_out[:, chan_ix] = pc_vecs[col_inds, chan_ix]
        test_mat = test_mat[col_inds]
        flips = np.sign(sim[row_inds, col_inds])
        flips[-1] *= np.prod(flips, axis=0)  # always flip an even number of factors
        for ix, f in enumerate(flips):
            align_out[ix, chan_ix] *= f
            test_mat[ix] *= f

        # Update ref_mat
        ref_mat = (ref_mat * (chan_ix + 1) + test_mat) / (chan_ix + 2)

    debug_plot(align_out)
    scipy.io.savemat(str(path)[:-4] + '_new.mat', {'pc_vecs': align_out})
