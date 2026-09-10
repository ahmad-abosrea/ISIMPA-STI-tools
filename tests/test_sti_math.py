# I-Simpa speech-intelligibility tools (STI, %ALcons).
# Copyright (C) 2026  Ahmad Abosrea
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
import math

import pytest

import sti_math as sm


def test_alpha_beta_are_self_consistent_at_perfect_transmission():
    """If every band's MTI = 1 (a perfect, distortion-free path), STI must
    be exactly 1. This only holds if alpha/beta are the correct matched
    pair -- it would fail immediately on a transcription typo in either
    table, so this is a real check, not a tautology."""
    mti = {f: 1.0 for f in sm.OCTAVE_BANDS_HZ}
    assert sm.compute_sti(mti) == pytest.approx(1.0, abs=1e-9)


def test_sti_is_zero_when_all_mti_zero():
    mti = {f: 0.0 for f in sm.OCTAVE_BANDS_HZ}
    assert sm.compute_sti(mti) == pytest.approx(0.0, abs=1e-9)


def test_sti_bounded_0_1_for_arbitrary_mti():
    mti = {125: 0.9, 250: 0.1, 500: 0.5, 1000: 0.3, 2000: 0.8, 4000: 0.2, 8000: 0.6}
    sti = sm.compute_sti(mti)
    assert 0.0 <= sti <= 1.0


def test_tcr_mtf_decreases_with_modulation_frequency():
    """Physical sanity: faster modulation is harder to preserve through a
    reverberant (low-pass-like) decay, so m(F) must decrease as F grows."""
    t60 = 1.0
    values = [sm.mtf_analytical_exponential_decay(t60, f) for f in sm.MODULATION_FREQS_HZ]
    assert all(values[i] > values[i + 1] for i in range(len(values) - 1))
    assert all(0.0 < v <= 1.0 for v in values)


def test_tcr_mtf_decreases_with_longer_t60():
    """Physical sanity: a longer reverberation time smears the signal more,
    so m(F) at a fixed modulation frequency must decrease as T60 grows."""
    f = 2.0
    m_short = sm.mtf_analytical_exponential_decay(0.4, f)
    m_long = sm.mtf_analytical_exponential_decay(2.5, f)
    assert m_long < m_short


def test_spps_echogram_mtf_matches_tcr_analytical_for_pure_exponential_decay():
    """The strongest available cross-check without a real GUI/echogram to
    test against: a purely exponential decay h^2(t) = exp(-t/tau) has a
    KNOWN closed-form MTF (the same formula the TCR branch uses). If the
    discrete echogram-based MTF integral is implemented correctly, running
    it on a finely-sampled synthetic exponential decay must reproduce the
    analytical TCR formula for the same tau, to within discretization
    error. This ties both branches together with an independent check that
    isn't just "the code doesn't crash"."""
    t60 = 1.2
    tau = t60 / 13.8
    dt = 0.0005  # 0.5 ms steps, similar order to a real SPPS time step
    duration = 8 * tau  # long enough that exp(-t/tau) has decayed to ~1e-3.5 or less
    n = int(duration / dt)
    times = [i * dt for i in range(1, n + 1)]
    h2 = [math.exp(-t / tau) for t in times]

    for f in sm.MODULATION_FREQS_HZ:
        m_echogram = sm.mtf_from_echogram(times, h2, f)
        m_analytical = sm.mtf_analytical_exponential_decay(t60, f)
        assert m_echogram == pytest.approx(m_analytical, rel=0.02), (
            f"mismatch at F={f} Hz: echogram={m_echogram}, analytical={m_analytical}"
        )


def test_apply_noise_snr_reduces_modulation_and_infinite_snr_is_a_no_op():
    m = 0.8
    assert sm.apply_noise_snr(m, snr_db=0.0) == pytest.approx(0.4, abs=1e-9)
    assert sm.apply_noise_snr(m, snr_db=100.0) == pytest.approx(m, abs=1e-6)


def test_transmission_index_bounds_and_monotonicity():
    assert sm.transmission_index(0.0) == pytest.approx(0.0, abs=1e-9)
    assert sm.transmission_index(1.0) == pytest.approx(1.0, abs=1e-9)
    assert sm.transmission_index(0.5) == pytest.approx(0.5, abs=1e-9)
    vals = [sm.transmission_index(m) for m in [0.1, 0.3, 0.5, 0.7, 0.9]]
    assert all(vals[i] < vals[i + 1] for i in range(len(vals) - 1))


def test_mti_for_band_averages_14_values():
    assert sm.mti_for_band([1.0] * 14) == pytest.approx(1.0)
    assert sm.mti_for_band([0.0] * 14) == pytest.approx(0.0)
    with pytest.raises(ValueError):
        sm.mti_for_band([1.0] * 13)


# --- Cross-validation against independent published references ---

ALPHA_MALE_IEC_TABLE_A1 = [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173]
BETA_MALE_IEC_TABLE_A1 = [0.085, 0.078, 0.065, 0.011, 0.047, 0.095]


def test_weighting_coefficients_match_published_iec_table_a1():
    """Independently confirmed against a published IEC 60268-16
    Edition 5.0 Table A.1 (male speech weighting + redundancy factors).
    Guards against a silent transcription error in the constants."""
    assert sm.ALPHA_MALE == ALPHA_MALE_IEC_TABLE_A1
    assert sm.BETA_MALE == BETA_MALE_IEC_TABLE_A1


def _becker_alcons_from_sti(sti):
    """Farrel Becker's widely used empirical STI -> %ALcons relation
    (used by DIRAC and others): %ALcons = 170.5405 * exp(-5.419 * STI)."""
    return 170.5405 * math.exp(-5.419 * sti)


@pytest.mark.parametrize("t60", [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
def test_sti_and_peutz_alcons_agree_via_independent_becker_relation(t60):
    """Cross-validates the two tools against each other using two
    INDEPENDENT literature sources:

      * this module's STI  (Houtgast/Steeneken MTF + IEC 60268-16 weighting)
        -> converted to %ALcons by Farrel Becker's empirical relation
      * Peutz's own saturated %ALcons asymptote (9 * T60), which is what
        alcons_tool falls back to beyond the critical distance rH

    These come from different derivations, so agreement is real evidence
    both are implemented correctly rather than a tautology. Measured
    agreement is within 22% relative across T60 = 1-4 s; the 25% threshold
    below is set just above that. Exact equality is NOT expected -- Becker's
    fit and Peutz's asymptote are two empirical approximations of the same
    underlying phenomenon. The threshold is still tight enough to catch a
    real implementation error (e.g. a factor-of-2 or a wrong tau)."""
    mti = {}
    for hz in sm.OCTAVE_BANDS_HZ:
        ti = [sm.transmission_index(sm.mtf_analytical_exponential_decay(t60, f))
              for f in sm.MODULATION_FREQS_HZ]
        mti[hz] = sm.mti_for_band(ti)
    sti = sm.compute_sti(mti)

    alcons_via_sti = _becker_alcons_from_sti(sti)
    alcons_peutz_saturated = 9.0 * t60

    rel_diff = abs(alcons_via_sti - alcons_peutz_saturated) / alcons_peutz_saturated
    assert rel_diff < 0.25, (
        "STI-derived %%ALcons (%.1f) and Peutz 9*T60 (%.1f) disagree by %.0f%% at T60=%.1fs"
        % (alcons_via_sti, alcons_peutz_saturated, rel_diff * 100, t60)
    )


def test_sti_decreases_monotonically_with_reverberation():
    """Sanity: more reverberation must never improve intelligibility."""
    values = []
    for t60 in [0.4, 0.8, 1.2, 2.0, 3.0, 4.0]:
        mti = {}
        for hz in sm.OCTAVE_BANDS_HZ:
            ti = [sm.transmission_index(sm.mtf_analytical_exponential_decay(t60, f))
                  for f in sm.MODULATION_FREQS_HZ]
            mti[hz] = sm.mti_for_band(ti)
        values.append(sm.compute_sti(mti))
    assert all(values[i] > values[i + 1] for i in range(len(values) - 1))


# --- Cross-confirmation against ANSI S3.5-1997 (SII) clause 5.2.3 ---
#
# ANSI S3.5-1997 defines the Speech Intelligibility Index (SII), NOT the STI.
# But its MTFI-based procedure (clause 5.2.3) shares the same core machinery
# as STI, so it is an INDEPENDENT standards-body confirmation of the three
# formulas below. Where they overlap, this module must match S3.5 exactly.
#
#   S3.5 eq (20):  M_f,i = |integral g(t) e^(-j2 pi f t) dt| / integral g(t) dt
#                  computed from the squared impulse response g(t)
#                  -> identical to mtf_from_echogram()
#   S3.5 eq (22):  R_f,i = 10 lg[ M_f,i / (1 - M_f,i) ]
#                  -> identical to modulation_to_apparent_snr_db()
#   S3.5 5.2.3.5:  R_f,i limited to the interval [-15, +15] dB
#                  -> identical to X_CLIP_DB
#
# S3.5 then DIVERGES: it averages R over NINE modulation frequencies
# (0.5, 1, 1.5, 2, 3, 4, 6, 8, 16 Hz) and feeds the result into the SII
# band-importance framework, whereas IEC 60268-16 uses FOURTEEN modulation
# frequencies (0.63-12.5 Hz) and the TI/MTI/alpha-beta framework. That is a
# difference between two different metrics, not an error in either.

ANSI_S35_MODULATION_FREQS_HZ = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 16.0]


def test_apparent_snr_formula_matches_ansi_s35_eq22():
    """R = 10 lg[m/(1-m)], clipped to [-15,+15] dB (ANSI S3.5 eq 22 and
    clause 5.2.3.5). Identical to the IEC 60268-16 X(F) step."""
    for m in [0.05, 0.2, 0.5, 0.8, 0.95]:
        expected = 10 * math.log10(m / (1 - m))
        expected = max(-15.0, min(15.0, expected))
        assert sm.modulation_to_apparent_snr_db(m) == pytest.approx(expected)
    # clipping at both rails
    assert sm.modulation_to_apparent_snr_db(0.999999) == pytest.approx(15.0)
    assert sm.modulation_to_apparent_snr_db(1e-9) == pytest.approx(-15.0)


def test_mtf_from_echogram_matches_ansi_s35_eq20_on_a_known_case():
    """S3.5 eq (20) is the same integral this module implements. Verified
    on a case with a closed-form answer: for g(t) = exp(-t/tau),
    |G(f)|/G(0) = 1/sqrt(1 + (2 pi f tau)^2)."""
    tau = 0.25
    dt = 0.0002
    times = [i * dt for i in range(1, int(10 * tau / dt) + 1)]
    g = [math.exp(-t / tau) for t in times]
    for f in ANSI_S35_MODULATION_FREQS_HZ:
        expected = 1.0 / math.sqrt(1 + (2 * math.pi * f * tau) ** 2)
        assert sm.mtf_from_echogram(times, g, f) == pytest.approx(expected, rel=0.02)


def test_iec_modulation_frequency_set_is_the_iec_one_not_the_ansi_one():
    """Guards the deliberate choice: this tool implements STI per
    IEC 60268-16 (14 modulation frequencies, 0.63-12.5 Hz), which is what
    was specified. ANSI S3.5's nine-frequency set belongs to SII and must
    not be substituted here."""
    assert len(sm.MODULATION_FREQS_HZ) == 14
    assert sm.MODULATION_FREQS_HZ[0] == 0.63
    assert sm.MODULATION_FREQS_HZ[-1] == 12.5
    assert sm.MODULATION_FREQS_HZ != ANSI_S35_MODULATION_FREQS_HZ


# --- rating scales (ISO 9921) ---

def test_iso9921_five_step_rating_boundaries():
    assert sm.iso9921_rating(0.80)[0] == "excellent"
    assert sm.iso9921_rating(0.75)[0] == "excellent"
    assert sm.iso9921_rating(0.70)[0] == "good"
    assert sm.iso9921_rating(0.60)[0] == "good"
    assert sm.iso9921_rating(0.50)[0] == "fair"
    assert sm.iso9921_rating(0.45)[0] == "fair"
    assert sm.iso9921_rating(0.35)[0] == "poor"
    assert sm.iso9921_rating(0.30)[0] == "poor"
    assert sm.iso9921_rating(0.29)[0] == "bad"
    assert sm.iso9921_rating(0.0)[0] == "bad"
    # the PB word score travels with the rating
    assert sm.iso9921_rating(0.80)[1] == "> 98 %"
    assert sm.iso9921_rating(0.0)[1] == "< 60 %"




def test_snr_is_speech_minus_noise_per_band():
    assert sm.snr_from_levels(65.0, 35.0) == pytest.approx(30.0)
    assert sm.snr_from_levels(40.0, 55.0) == pytest.approx(-15.0)


def test_noise_lowers_sti_and_noise_free_is_the_upper_bound():
    """Adding background noise can only reduce STI, so the noise-free
    result is an upper bound -- which is what the tool reports when no
    noise spectrum is given."""
    mti_free, mti_noisy = {}, {}
    for hz in sm.OCTAVE_BANDS_HZ:
        ti_free, ti_noisy = [], []
        for f in sm.MODULATION_FREQS_HZ:
            m = sm.mtf_analytical_exponential_decay(1.0, f)
            ti_free.append(sm.transmission_index(sm.apply_noise_snr(m, 99.0)))
            ti_noisy.append(sm.transmission_index(sm.apply_noise_snr(m, 5.0)))
        mti_free[hz] = sm.mti_for_band(ti_free)
        mti_noisy[hz] = sm.mti_for_band(ti_noisy)
    assert sm.compute_sti(mti_noisy) < sm.compute_sti(mti_free)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


def test_letter_qualification_bands_stay_removed():
    """The A+..U letter scale was removed because its thresholds could not
    be verified against IEC 60268-16 (paywalled), and neither the a published IEC 60268-16 guide
    IEC 60268-16 guide nor a published Edition 5 compliance reference uses
    letter bands. This test stops it being reintroduced by habit."""
    assert not hasattr(sm, "IEC_QUALIFICATION_BANDS")
    assert not hasattr(sm, "iec_qualification_band")


def test_five_step_boundaries_match_the_published_guide():
    """0.30 / 0.45 / 0.60 / 0.75 -- independently confirmed against the
    published IEC 60268-16 guide's quality-category table."""
    assert sm.iso9921_rating(0.20)[0] == "bad"
    assert sm.iso9921_rating(0.30)[0] == "poor"
    assert sm.iso9921_rating(0.44)[0] == "poor"
    assert sm.iso9921_rating(0.45)[0] == "fair"
    assert sm.iso9921_rating(0.59)[0] == "fair"
    assert sm.iso9921_rating(0.60)[0] == "good"
    assert sm.iso9921_rating(0.74)[0] == "good"
    assert sm.iso9921_rating(0.75)[0] == "excellent"


def test_alpha_beta_match_two_independent_published_sources():
    """Two independent published references to IEC 60268-16
    Edition 5 compliance" page list exactly these male alpha/beta values."""
    assert sm.ALPHA_MALE == [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173]
    assert sm.BETA_MALE == [0.085, 0.078, 0.065, 0.011, 0.047, 0.095]
    assert sm.OCTAVE_BANDS_HZ == [125, 250, 500, 1000, 2000, 4000, 8000]
    assert len(sm.MODULATION_FREQS_HZ) == 14


def test_modulation_frequencies_are_the_published_nominal_values():
    """The 14 nominal modulation frequencies, identical in the a published IEC 60268-16 guide
    IEC 60268-16 guide and a published Edition 5 compliance reference. A 6.25
    typo for 6.30 lived here until that reference's list caught it."""
    assert sm.MODULATION_FREQS_HZ == [
        0.63, 0.80, 1.00, 1.25, 1.60, 2.00, 2.50, 3.15,
        4.00, 5.00, 6.30, 8.00, 10.00, 12.50]


def test_nominal_frequencies_track_the_exact_generator():
    """They are generated as 10^((-2:11)/10); the published nominal
    values are those rounded, so each must sit within ~1% of exact."""
    for got, k in zip(sm.MODULATION_FREQS_HZ, range(-2, 12)):
        exact = 10 ** (k / 10.0)
        assert abs(got - exact) / exact < 0.012, (got, exact)


# ---------- auditory masking + absolute reception threshold ----------
#
# These constants cannot be checked against the paywalled standard, so the
# tests pin them against TWO independent published sources that agree:
#   [1] Steeneken, "Basics of the STI measuring method", Tables I and II
#   [2] Architectural Institute of Japan benchmark A-21, which reproduces
#       IEC 60268-16 Table A.2 (slope of masking) and Table A.3 (reception
#       threshold) verbatim, attributed to IEC 60268-16:2020

def test_threshold_values_match_both_sources():
    """Steeneken Table II and IEC Table A.3 are identical."""
    assert sm.ABSOLUTE_RECEPTION_THRESHOLD_DB == [46.0, 27.0, 12.0, 6.5, 7.5, 8.0, 12.0]
    assert len(sm.ABSOLUTE_RECEPTION_THRESHOLD_DB) == len(sm.OCTAVE_BANDS_HZ)


def test_masking_slope_matches_iec_table_a2():
    """The four segments of the piecewise slope of masking specified in
    IEC 60268-16 Annex A, Table A.2, as reproduced in the Architectural
    Institute of Japan benchmark A-21."""
    assert sm.masking_slope_db(50.0) == pytest.approx(0.5 * 50.0 - 65.0)
    assert sm.masking_slope_db(65.0) == pytest.approx(1.8 * 65.0 - 146.9)
    assert sm.masking_slope_db(80.0) == pytest.approx(0.5 * 80.0 - 59.8)
    assert sm.masking_slope_db(120.0) == -10.0


def test_masking_slope_is_continuous_at_both_joins():
    """A transcription slip in any coefficient would open a step here."""
    for edge in (63.0, 67.0):
        below = sm.masking_slope_db(edge - 1e-6)
        at = sm.masking_slope_db(edge)
        assert abs(below - at) < 1e-3, (edge, below, at)


@pytest.mark.parametrize("lo,hi,nominal_slope,published_amf", [
    (46, 55, -40, 0.000100),
    (56, 65, -35, 0.000316),
    (66, 75, -25, 0.003162),
    (76, 85, -20, 0.010000),
    (86, 95, -15, 0.031622),
])
def test_continuous_formula_agrees_with_steenekens_banded_table(
        lo, hi, nominal_slope, published_amf):
    """Source [2] is the standard's continuous ramp; source [1] is the
    banded empirical data it was smoothed from. Each band spans about
    bands. Across each band the continuous value must straddle the
    published nominal slope."""
    at_lo = sm.masking_slope_db(float(lo))
    at_hi = sm.masking_slope_db(float(hi))
    assert min(at_lo, at_hi) <= nominal_slope <= max(at_lo, at_hi), (at_lo, at_hi)
    assert 10 ** (nominal_slope / 10.0) == pytest.approx(published_amf, rel=0.001)


def test_amf_is_ten_to_the_slope_over_ten():
    """Steeneken's Table I publishes both; they must be consistent."""
    for level in (50.0, 60.0, 70.0, 80.0, 90.0, 110.0):
        assert sm.auditory_masking_factor(level) == pytest.approx(
            10 ** (sm.masking_slope_db(level) / 10.0))


def test_masking_off_reproduces_the_plain_noise_term_exactly():
    """Switching the hearing model off must leave the previous behaviour
    bit-for-bit, so the option cannot silently change old results."""
    speech = {hz: 60.0 for hz in sm.OCTAVE_BANDS_HZ}
    noise = {hz: 40.0 for hz in sm.OCTAVE_BANDS_HZ}
    factors = sm.modulation_reduction_factors(speech, noise, use_masking=False)
    for hz in sm.OCTAVE_BANDS_HZ:
        assert factors[hz] == pytest.approx(sm.apply_noise_snr(1.0, 20.0), abs=1e-12)


def test_masking_can_only_reduce_the_modulation_index():
    speech = {hz: 70.0 for hz in sm.OCTAVE_BANDS_HZ}
    off = sm.modulation_reduction_factors(speech, None, use_masking=False)
    on = sm.modulation_reduction_factors(speech, None, use_masking=True)
    for hz in sm.OCTAVE_BANDS_HZ:
        assert on[hz] <= off[hz] + 1e-12
        assert 0.0 < on[hz] <= 1.0


def test_lowest_band_is_never_masked_only_thresholded():
    """Masking comes from the band BELOW, and 125 Hz has none. Making the
    250 Hz band deafening must not change the 125 Hz factor at all."""
    quiet = {hz: 60.0 for hz in sm.OCTAVE_BANDS_HZ}
    loud_above = dict(quiet)
    loud_above[250] = 110.0
    a = sm.modulation_reduction_factors(quiet, None, True)
    b = sm.modulation_reduction_factors(loud_above, None, True)
    assert a[125] == pytest.approx(b[125])
    assert b[500] < a[500]          # 250 Hz now masks 500 Hz


def test_a_loud_low_band_masks_the_band_above_it():
    base = {hz: 55.0 for hz in sm.OCTAVE_BANDS_HZ}
    boomy = dict(base)
    boomy[125] = 95.0
    assert (sm.modulation_reduction_factors(boomy, None, True)[250]
            < sm.modulation_reduction_factors(base, None, True)[250])


def test_threshold_dominates_when_speech_is_near_inaudible():
    """At 125 Hz the threshold is 46 dB, so a 30 dB signal must be crushed
    while a 90 dB one is barely touched."""
    faint = sm.modulation_reduction_factors({hz: 30.0 for hz in sm.OCTAVE_BANDS_HZ},
                                            None, True)
    loud = sm.modulation_reduction_factors({hz: 90.0 for hz in sm.OCTAVE_BANDS_HZ},
                                           None, True)
    assert faint[125] < 0.05
    assert loud[125] > 0.99
