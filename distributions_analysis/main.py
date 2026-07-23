import csv
import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.signal import savgol_filter

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BEFORE_CSV = os.path.join(SCRIPT_DIR, "distribution_before.csv")
AFTER_CSV = os.path.join(SCRIPT_DIR, "distribution_after.csv")
SHEARED_CSV = os.path.join(SCRIPT_DIR, "distribution_sheared.csv")
SPRI_CSV = os.path.join(SCRIPT_DIR, "distribution_spri.csv")
SPRI_SHEARED_RATIO_CSV = os.path.join(SCRIPT_DIR, "spri_sheared_ratio.csv")
PROBABILITY_OF_KEEPING_CSV = os.path.join(SCRIPT_DIR, "probability_of_keeping.csv")
MODIFIED_BEFORE_CSV = os.path.join(SCRIPT_DIR, "distribution_before_modified.csv")
OLD_FRAGMENT_SIZE_DIST_TXT = os.path.join(SCRIPT_DIR, "..", "radiSeqData", "FragmentSizeDist_induceSeq_old.txt")
NEW_FRAGMENT_SIZE_DIST_TXT = os.path.join(SCRIPT_DIR, "..", "radiSeqData", "FragmentSizeDist_induceSeq.txt")


def load_distribution(csv_path):
    """Loads a distribution csv file (with a 'bp,FU' header) into a pair of numpy arrays (bp, FU)."""
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
    bp = data[:, 0]
    fu = data[:, 1]
    return bp, fu


def load_fragment_size_dist_txt(txt_path):
    """
    Loads a whitespace-separated 'length count' fragment size distribution file with no header
    (the format used by the C++ readFragmentSizeDist function, e.g. FragmentSizeDist_induceSeq_old.txt),
    into a pair of numpy arrays (bp, count).
    """
    data = np.loadtxt(txt_path)
    bp = data[:, 0]
    count = data[:, 1]
    return bp, count


def graph_distributions(normalize=False):
    """
    Loads the before/after distributions and plots them together against fragment size.
    If normalize is True, each distribution's FU values are divided by their own sum first,
    so that each one sums to 1 (comparable by shape rather than absolute intensity).
    """
    before_bp, before_fu = load_distribution(BEFORE_CSV)
    after_bp, after_fu = load_distribution(AFTER_CSV)

    if normalize:
        before_fu = before_fu / before_fu.sum()
        after_fu = after_fu / after_fu.sum()

    plt.plot(before_bp, before_fu, color="green", label="before")
    plt.plot(after_bp, after_fu, color="orange", label="after")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Normalized proportion" if normalize else "Sample Intensity [Normalized FU]")
    plt.title("Fragment size distributions")
    plt.legend()
    plt.grid(True)
    plt.show()

def graph_shifted(shift, normalize=False, scale_after=1):
    """
    Loads the before/after distributions and plots them together, with the before distribution's
    bp values shifted by 'shift'. If normalize is True, each distribution's FU values are divided
    by their own sum first, so that each one sums to 1 (comparable by shape rather than absolute
    intensity). scale_after multiplies the after distribution's FU values (applied after
    normalization, if normalize is True).
    """
    before_bp, before_fu = load_distribution(BEFORE_CSV)
    after_bp, after_fu = load_distribution(AFTER_CSV)

    if normalize:
        before_fu = before_fu / before_fu.sum()
        after_fu = after_fu / after_fu.sum()

    after_fu = after_fu * scale_after

    plt.plot(before_bp + shift, before_fu, color="green", label=f"before (shifted by {shift})")
    plt.plot(after_bp, after_fu, color="orange", label=f"after (scaled by {scale_after})")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Normalized proportion" if normalize else "Sample Intensity [Normalized FU]")
    plt.title("Fragment size distributions (before shifted)")
    plt.legend()
    plt.grid(True)
    plt.show()


def graph_sheared_spri(scale_spri=1):
    """
    Loads the sheared/spri distributions and plots them together against fragment size.
    scale_spri multiplies the spri distribution's FU values.
    """
    sheared_bp, sheared_fu = load_distribution(SHEARED_CSV)
    spri_bp, spri_fu = load_distribution(SPRI_CSV)

    spri_fu = spri_fu * scale_spri

    plt.plot(sheared_bp, sheared_fu, color="red", label="sheared")
    plt.plot(spri_bp, spri_fu, color="navy", label=f"spri (scaled by {scale_spri})")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Sample Intensity [Normalized FU]")
    plt.title("Sheared vs SPRI (1.2x) fragment size distributions")
    plt.legend()
    plt.grid(True)
    plt.show()


def graph_spri_sheared_ratio(scale=1, window=21):
    """
    Loads the sheared/spri distributions, scales the spri distribution's FU values by 'scale',
    then divides the (scaled) spri distribution by the sheared distribution at each of the spri
    distribution's bp values (the sheared distribution is linearly interpolated onto that bp
    grid, since the two were extracted independently and don't share the same bp points).

    Smooths the resulting quotient with a Savitzky-Golay filter, then holds the smoothed
    quotient at 1 for the rest of the range once it first reaches 1, and extends the range
    down to 0 bp with a ratio of 0 there (both distributions are essentially zero below the
    lowest extracted bp value). Saves the smoothed quotient (paired with bp) to
    SPRI_SHEARED_RATIO_CSV, and graphs the smoothed and unsmoothed quotient together against
    fragment size.
    """
    sheared_bp, sheared_fu = load_distribution(SHEARED_CSV)
    spri_bp, spri_fu = load_distribution(SPRI_CSV)

    spri_fu = spri_fu * scale
    sheared_fu_interp = np.interp(spri_bp, sheared_bp, sheared_fu)
    ratio = spri_fu / sheared_fu_interp

    # Extend the range down to 0 bp, with the ratio being 0 there, before smoothing: this lets
    # the filter blend gradually into the zero baseline instead of smoothing the two pieces
    # separately and stitching them together (which left a jump/dip artifact at the seam).
    bp_step = spri_bp[1] - spri_bp[0]
    zero_bp = np.arange(0, spri_bp[0], bp_step)
    spri_bp = np.concatenate([zero_bp, spri_bp])
    ratio = np.concatenate([np.zeros_like(zero_bp), ratio])

    # Savitzky-Golay window must be odd and no longer than the data itself.
    window_length = min(window, len(ratio) - (1 - len(ratio) % 2))
    ratio_smoothed = savgol_filter(ratio, window_length=window_length, polyorder=3)
    ratio_smoothed = np.clip(ratio_smoothed, 0, None)

    # Once the smoothed ratio first reaches 1, hold it at 1 for the rest of the range.
    reaches_one = np.where(ratio_smoothed >= 1)[0]
    if len(reaches_one) > 0:
        ratio_smoothed[reaches_one[0]:] = 1

    with open(SPRI_SHEARED_RATIO_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["bp", "quotient"])
        for b, q in zip(spri_bp, ratio_smoothed):
            writer.writerow([round(float(b), 1), round(float(q), 4)])

    plt.plot(spri_bp, ratio, color="plum", label="unsmoothed")
    plt.plot(spri_bp, ratio_smoothed, color="purple", label="smoothed")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Quotient")
    plt.title(f"SPRI (1.2x, scaled by {scale}) / Sheared ratio")
    plt.legend()
    plt.grid(True)
    plt.show()


def compute_probability_of_keeping():
    """
    Computes, for every integer fragment size x from 1 to 2000 bp, the probability of keeping
    a fragment of that size: r(x + 25 + 58) * r(x - 25 + 58 + 6) * 0.97**2, where r(L) is the
    spri_sheared_ratio value at L (linearly interpolated from SPRI_SHEARED_RATIO_CSV). Saves the
    resulting (bp, probability) pairs to PROBABILITY_OF_KEEPING_CSV.
    """
    ratio_bp, ratio_q = load_distribution(SPRI_SHEARED_RATIO_CSV)

    def r(L):
        return np.interp(L, ratio_bp, ratio_q)

    x = np.arange(1, 2001)
    probability = r(x + 25 + 58) * r(x - 25 + 58 + 6) * 0.97**2

    with open(PROBABILITY_OF_KEEPING_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["bp", "probability"])
        for b, p in zip(x, probability):
            writer.writerow([int(b), round(float(p), 6)])


def graph_probability_of_keeping():
    """Loads and plots the probability_of_keeping distribution against fragment size."""
    bp, probability = load_distribution(PROBABILITY_OF_KEEPING_CSV)

    plt.plot(bp, probability, color="teal")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Probability of keeping")
    plt.title("Probability of keeping vs fragment size")
    plt.grid(True)
    plt.show()


def compute_modified_before(low_bp=120, high_bp=1000, taper_width=30, concave_fit_start=570, concave_fit_end=822):
    """
    Loads distribution_before and creates a modified version of it that is 0 below low_bp,
    smoothly (rather than abruptly) reaches 0 as bp approaches low_bp from above, is unchanged
    in the middle of the range, and smoothly reaches 0 again as bp approaches high_bp from below,
    staying 0 above high_bp (removing the small flat tail distribution_before has above high_bp).

    Before doing that, it also smooths out a small bump the raw tail has around 950bp (where the
    distribution briefly rises/plateaus instead of continuing to smoothly decay): an exponential
    decay A*exp(-k*bp) is fit (in log space, via linear regression) to the clean, smoothly
    decreasing region [concave_fit_start, concave_fit_end], and used to replace all values above
    concave_fit_end, so the curve keeps decaying the way it already does from concave_fit_start
    onwards, all the way down to where the high_bp taper brings it to 0.

    This is done by multiplying the FU values by a window that is 0 outside [low_bp, high_bp] and
    ramps smoothly (via a smoothstep) from 0 to 1 over the first taper_width bp above low_bp, and
    from 1 to 0 over the last taper_width bp below high_bp.

    Saves the result to MODIFIED_BEFORE_CSV and returns (bp, modified_fu).
    """
    before_bp, before_fu = load_distribution(BEFORE_CSV)

    # Replace the noisy tail above concave_fit_end with the smooth exponential decay fit to the
    # clean region just before it, removing the small bump/plateau it otherwise has around 950bp.
    fit_region = (before_bp >= concave_fit_start) & (before_bp <= concave_fit_end)
    decay_rate, log_amplitude = np.polyfit(before_bp[fit_region], np.log(before_fu[fit_region]), 1)
    replace_region = before_bp > concave_fit_end
    before_fu[replace_region] = np.exp(log_amplitude) * np.exp(decay_rate * before_bp[replace_region])

    def smoothstep(t):
        t = np.clip(t, 0, 1)
        return 3 * t**2 - 2 * t**3

    rising = np.where(before_bp <= low_bp, 0.0, smoothstep((before_bp - low_bp) / taper_width))
    falling = np.where(before_bp >= high_bp, 0.0, smoothstep((high_bp - before_bp) / taper_width))
    window = rising * falling

    modified_fu = before_fu * window

    with open(MODIFIED_BEFORE_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["bp", "FU"])
        for b, f_val in zip(before_bp, modified_fu):
            writer.writerow([round(float(b), 1), round(float(f_val), 4)])

    return before_bp, modified_fu


def graph_modified_before():
    """Loads and plots the modified before distribution, together with the original, against fragment size."""
    before_bp, before_fu = load_distribution(BEFORE_CSV)
    bp, fu = load_distribution(MODIFIED_BEFORE_CSV)

    plt.plot(before_bp, before_fu, color="limegreen", label="before")
    plt.plot(bp, fu, color="darkgreen", label="modified before")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Sample Intensity [Normalized FU]")
    plt.title("Modified before distribution")
    plt.legend()
    plt.grid(True)
    plt.show()


def graph_modified_before_vs_old_dist():
    """
    Loads the old induce_seq fragment size distribution (FragmentSizeDist_induceSeq_old.txt, a
    whitespace-separated 'length count' file with no header, unlike the CSV distributions used
    elsewhere in this file) and the modified before distribution, normalizes each to sum to 1,
    and plots them together against fragment size.
    """
    old_bp, old_count = load_fragment_size_dist_txt(OLD_FRAGMENT_SIZE_DIST_TXT)
    modified_bp, modified_fu = load_distribution(MODIFIED_BEFORE_CSV)

    old_count = old_count / old_count.sum()
    modified_fu = modified_fu / modified_fu.sum()

    plt.plot(old_bp, old_count, color="steelblue", label="old induce_seq fragment size dist")
    plt.plot(modified_bp, modified_fu, color="darkgreen", label="modified before")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Normalized proportion")
    plt.title("Modified before vs old induce_seq fragment size distribution")
    plt.legend()
    plt.grid(True)
    plt.show()


def save_modified_before_as_fragment_size_dist(low_bp=120, high_bp=1000):
    """
    Loads the modified before distribution and writes it out in the same whitespace-separated
    'length count' format as FragmentSizeDist_induceSeq_old.txt (see load_fragment_size_dist_txt),
    to NEW_FRAGMENT_SIZE_DIST_TXT, so it can be used as the induce_seq fragment size distribution
    file. One row is written for every integer bp from 1 to 2000, with the count linearly
    interpolated from the modified before distribution. low_bp/high_bp should match the values
    used in compute_modified_before: they are re-applied here to force exact 0 outside that range,
    since interpolating onto the integer grid can otherwise leak tiny nonzero values across the
    large gap the source distribution has just below low_bp (where its marker peak was excluded).
    """
    bp, fu = load_distribution(MODIFIED_BEFORE_CSV)

    lengths = np.arange(1, 2001)
    counts = np.interp(lengths, bp, fu)
    counts[(lengths <= low_bp) | (lengths >= high_bp)] = 0.0

    with open(NEW_FRAGMENT_SIZE_DIST_TXT, "w") as f:
        for length, count in zip(lengths, counts):
            f.write(f"{int(length)} {round(float(count), 4)}\n")


def _modified_before(before_bp, before_fu, x, shift, percent_shifted, percent_half_shifted):
    """
    Evaluates, at bp values x, the modified-before mixture: percent_shifted of the before
    distribution shifted by shift, percent_half_shifted of it shifted by shift/2, and the
    remainder (1 - percent_shifted - percent_half_shifted) left unshifted.
    """
    def shifted_before(amount):
        # Interpolate the before distribution as if its bp values were all moved by 'amount';
        # 0 outside the shifted distribution's domain, since there is no data to extrapolate from there.
        return np.interp(x - amount, before_bp, before_fu, left=0, right=0)

    percent_unshifted = 1 - percent_shifted - percent_half_shifted
    return (
        percent_unshifted * shifted_before(0)
        + percent_half_shifted * shifted_before(shift / 2)
        + percent_shifted * shifted_before(shift)
    )


def _plot_fit(after_bp, model_before_fu, scaled_after_fu, shift, percent_shifted, percent_half_shifted, scale):
    plt.plot(after_bp, model_before_fu, color="green", label="modified before")
    plt.plot(after_bp, scaled_after_fu, color="orange", label="scaled after")
    plt.xscale("log")
    plt.xlabel("Size [bp]")
    plt.ylabel("Sample Intensity [Normalized FU]")
    plt.title(
        f"shift={shift}, percent_shifted={percent_shifted:.3f}, "
        f"percent_half_shifted={percent_half_shifted:.3f}, scale={scale:.3f}"
    )
    plt.legend()
    plt.grid(True)
    plt.show()


def fit_parameters(shift, initial_guess=(0.9, 0.05, 1.0)):
    """
    Finds the (percent_shifted, percent_half_shifted, scale) that minimize the least-squares
    distance between a scaled after distribution and a modified before distribution (see
    _modified_before), subject to percent_shifted > 0.8. shift is a fixed parameter of the
    function (in bp), not something being fit. initial_guess is the starting
    (percent_shifted, percent_half_shifted, scale) point for the optimizer.

    Graphs the fitted modified-before distribution against the fitted scaled-after distribution,
    with shift and the three fitted parameter values shown on the plot. Returns
    (percent_shifted, percent_half_shifted, scale).
    """
    before_bp, before_fu = load_distribution(BEFORE_CSV)
    after_bp, after_fu = load_distribution(AFTER_CSV)

    # The least-squares fit is only calculated over the 400-1000bp region (the graph afterwards
    # still shows the whole range, using the fitted parameters).
    fit_region = (after_bp >= 410) & (after_bp <= 720)
    after_bp_fit, after_fu_fit = after_bp[fit_region], after_fu[fit_region]

    def objective(params):
        percent_shifted, percent_half_shifted, scale = params
        # Keep the two percentages a valid split of the before distribution's mass (each in
        # [0,1], summing to at most 1), percent_shifted above 0.8, and the scale non-negative;
        # penalize the rest.
        if (
            percent_shifted <= 0.6
            or percent_half_shifted < 0
            or percent_shifted + percent_half_shifted > 1
            or scale < 0
        ):
            return 1e18
        model = _modified_before(before_bp, before_fu, after_bp_fit, shift, percent_shifted, percent_half_shifted)
        return np.sum((scale * after_fu_fit - model) ** 2)

    # Nelder-Mead (derivative-free) is used instead of a gradient-based method: the objective's
    # gradient w.r.t. scale is orders of magnitude larger than w.r.t. the two percentages (since
    # FU values are in the hundreds), which made gradient-based solvers (e.g. SLSQP) stall at
    # the initial guess instead of taking a step.
    result = minimize(objective, initial_guess, method="Nelder-Mead")
    percent_shifted, percent_half_shifted, scale = result.x

    model_before_fu = _modified_before(before_bp, before_fu, after_bp, shift, percent_shifted, percent_half_shifted)
    _plot_fit(after_bp, model_before_fu, scale * after_fu, shift, percent_shifted, percent_half_shifted, scale)

    return percent_shifted, percent_half_shifted, scale


def graph_with_parameters(shift, percent_shifted, percent_half_shifted, scale):
    """
    Graphs the modified-before distribution and the scaled-after distribution (see
    _modified_before) using manually-specified parameter values, instead of fitting them.
    """
    before_bp, before_fu = load_distribution(BEFORE_CSV)
    after_bp, after_fu = load_distribution(AFTER_CSV)

    model_before_fu = _modified_before(before_bp, before_fu, after_bp, shift, percent_shifted, percent_half_shifted)
    _plot_fit(after_bp, model_before_fu, scale * after_fu, shift, percent_shifted, percent_half_shifted, scale)



if __name__ == "__main__":
    # graph_distributions(normalize=True)
    # graph_modified_before()
    graph_modified_before_vs_old_dist()
    # graph_shifted(58 * 2 + 6, normalize=True, scale_after=0.77)
    # graph_distributions(normalize=True)
    # graph_shifted(58*2+6, scale_after=0.95)
    # graph_sheared_spri(scale_spri=0.67)
    # graph_spri_sheared_ratio(scale=0.67, window=101)
    # graph_probability_of_keeping()
    # fit_parameters(58*2+6, initial_guess=(0.9, 0.05, 1))
    # graph_with_parameters(58*2+6, 0.5, 0.1, 0.95)
