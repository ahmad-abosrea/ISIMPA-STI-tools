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
"""
Pure, dependency-free IEC 60268-16 (full STI, male-voice weighting) math.
Deliberately separated from the I-Simpa scripting glue (sti_tool/__init__.py,
which must live inside I-Simpa-main/src/isimpa/resources/UserScript/ and
can't be imported/pytested standalone) so this can be validated on its own
outside the GUI.

Constants below (IEC 60268-16, male talker octave weighting alpha_i and
adjacent-band redundancy beta_i) are the standard published
Houtgast/Steeneken coefficients. Self-consistency check baked into
test_sti_math.py: when every band's MTI = 1 (a perfect, distortion-free
transmission path), STI must equal exactly 1 -- this only holds if alpha
and beta are the correct, matched pair, so it is a strong check against a
transcription error in these numbers, not just a sanity bound.
"""
from __future__ import annotations

import cmath
import math
from typing import Dict, List, Sequence

OCTAVE_BANDS_HZ: List[int] = [125, 250, 500, 1000, 2000, 4000, 8000]
MODULATION_FREQS_HZ: List[float] = [
    0.63, 0.80, 1.00, 1.25, 1.60, 2.00, 2.50, 3.15,
    4.00, 5.00, 6.30, 8.00, 10.00, 12.50,
]

# Male-voice weighting (alpha, one per octave band) and adjacent-band
# redundancy correction (beta, one per adjacent pair -- 6 for 7 bands).
ALPHA_MALE: List[float] = [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173]
BETA_MALE: List[float] = [0.085, 0.078, 0.065, 0.011, 0.047, 0.095]

X_CLIP_DB = 15.0


def mtf_analytical_exponential_decay(t60_s: float, f_mod_hz: float) -> float:
    """TCR branch: statistical model, no real impulse response. Analytical
    MTF for a purely exponential reverberant decay.
        tau = T60 / 13.8   (13.8 = 60/ln(10)*ln(10) ... standard T60->tau constant, i.e. 60/(10*log10(e))... )
        m(F) = 1 / sqrt(1 + (2*pi*F*tau)^2)
    """
    if t60_s <= 0:
        raise ValueError("T60 must be positive")
    tau = t60_s / 13.8
    return 1.0 / math.sqrt(1.0 + (2 * math.pi * f_mod_hz * tau) ** 2)


def mtf_from_echogram(times_s: Sequence[float], h_squared: Sequence[float], f_mod_hz: float) -> float:
    """SPPS branch: true echogram available. Discrete numerical version of
        m(F) = | integral( h^2(t) * exp(-j*2*pi*F*t) dt ) | / integral( h^2(t) dt )
    Samples need not be uniformly spaced -- each sample's local time-step
    width is used as its integration weight, which is why weighting is
    computed per-sample rather than assumed constant. This also means an
    overall constant dt cancels between numerator and denominator, but is
    still computed explicitly here for correctness with non-uniform
    sampling and for auditability.
    """
    n = len(times_s)
    if n != len(h_squared):
        raise ValueError("times_s and h_squared must be the same length")
    if n < 2:
        raise ValueError("Need at least 2 echogram samples")

    weights = _sample_weights(times_s)

    numerator = 0j
    denominator = 0.0
    for t, h2, w in zip(times_s, h_squared, weights):
        numerator += h2 * w * cmath.exp(-1j * 2 * math.pi * f_mod_hz * t)
        denominator += h2 * w
    if denominator <= 0:
        raise ValueError("Echogram has zero total energy")
    return abs(numerator) / denominator


def _sample_weights(times_s: Sequence[float]) -> List[float]:
    """Standard trapezoidal-rule weight per sample (works for uniform or
    non-uniform sampling): interior samples get (t[i+1]-t[i-1])/2,
    endpoints get half of their single adjacent interval."""
    n = len(times_s)
    weights = [0.0] * n
    weights[0] = (times_s[1] - times_s[0]) / 2.0
    weights[-1] = (times_s[-1] - times_s[-2]) / 2.0
    for i in range(1, n - 1):
        weights[i] = (times_s[i + 1] - times_s[i - 1]) / 2.0
    return weights


def apply_noise_snr(m: float, snr_db: float) -> float:
    """m'(F) = m(F) / (1 + 10^(-SNR/10))"""
    return m / (1.0 + 10 ** (-snr_db / 10.0))


def modulation_to_apparent_snr_db(m_prime: float) -> float:
    """X(F) = 10*log10( m'(F) / (1 - m'(F)) ), clipped to [-15, +15] dB"""
    if m_prime >= 1.0:
        return X_CLIP_DB
    if m_prime <= 0.0:
        return -X_CLIP_DB
    x = 10 * math.log10(m_prime / (1.0 - m_prime))
    return max(-X_CLIP_DB, min(X_CLIP_DB, x))


def transmission_index(m_prime: float) -> float:
    """TI(F) = (X(F) + 15) / 30, in [0, 1]"""
    x = modulation_to_apparent_snr_db(m_prime)
    return (x + X_CLIP_DB) / (2 * X_CLIP_DB)


def mti_for_band(ti_values_14: Sequence[float]) -> float:
    """MTI(band) = average of TI(F) over the 14 modulation frequencies."""
    if len(ti_values_14) != 14:
        raise ValueError("Expected exactly 14 TI values (one per modulation frequency)")
    return sum(ti_values_14) / 14.0


def compute_sti(
    mti_by_octave_hz: Dict[int, float],
    alpha: Sequence[float] = ALPHA_MALE,
    beta: Sequence[float] = BETA_MALE,
) -> float:
    """STI = sum(alpha_i * MTI_i) - sum(beta_j * sqrt(MTI_j * MTI_j+1))
    over the 7 octave bands (in OCTAVE_BANDS_HZ order) and their 6 adjacent
    pairs. Clipped to [0, 1] for numerical safety."""
    mti_ordered = [mti_by_octave_hz[f] for f in OCTAVE_BANDS_HZ]
    if len(mti_ordered) != 7:
        raise ValueError("Expected MTI for all 7 IEC 60268-16 octave bands")

    weighted_sum = sum(a * m for a, m in zip(alpha, mti_ordered))
    redundancy = sum(
        b * math.sqrt(max(0.0, mti_ordered[i]) * max(0.0, mti_ordered[i + 1]))
        for i, b in enumerate(beta)
    )
    sti = weighted_sum - redundancy
    return max(0.0, min(1.0, sti))


# ---------------------------------------------------------------------------
# Intelligibility rating scales
# ---------------------------------------------------------------------------
# Two published scales, both reported because they answer different questions.
#
# 1. ISO 9921 five-step rating, with the associated
#    phonetically-balanced word score.
# The lettered "qualification band" scale (A+ .. U) that used to be here has
# been REMOVED. It came from a secondary source (a vendor application note
# citing IEC 60268-16 Annex F) and its thresholds could not be verified
# against the normative text, which is paywalled. Neither the a published IEC 60268-16 guide
# IEC 60268-16 guide nor a published Edition 5 compliance reference
# reference uses letter bands. Reporting one labelled "IEC Annex F" claimed
# more authority than could be backed, so only the five-step rating is kept.

ISO9921_RATINGS = [          # (lower bound, label, PB word score)
    (0.75, "excellent", "> 98 %"),
    (0.60, "good", "93 to 98 %"),
    (0.45, "fair", "80 to 93 %"),
    (0.30, "poor", "60 to 80 %"),
    (0.00, "bad", "< 60 %"),
]

STI_JND = 0.03  # Bradley et al.


def iso9921_rating(sti):
    """ISO 9921 five-step intelligibility rating. Returns (label, pb_score)."""
    for lower, label, pb in ISO9921_RATINGS:
        if sti >= lower:
            return label, pb
    return "bad", "< 60 %"


# iec_qualification_band() was removed along with the letter-band table.


def snr_from_levels(speech_spl_db, noise_spl_db):
    """SNR per octave band = speech level - background noise level, both as
    SPL at the receiver ('Signal-to-noise ratio -- difference
    between the sound pressure level of the speech or test signal and the
    sound pressure level of the background noise')."""
    return speech_spl_db - noise_spl_db


# ---------------------------------------------------------------------------
# Noise Criterion (NC) background-noise spectra
# ---------------------------------------------------------------------------
# Octave-band SPL limits in dB at 63 .. 8000 Hz. Beranek's NC curves, as
# tabulated in the ASHRAE Handbook (Ch. 48) / ANSI-ASA S12.2. Values were
# cross-checked against two independent published tables which agreed
# exactly across NC-15 .. NC-50; the NC-55 .. NC-65 rows come from the
# source that tabulates the full range.
#
# The NC rating of a room is the highest curve "touched" by its measured
# spectrum, so using a curve as a background-noise input represents the
# worst-case spectrum still qualifying for that rating.
#
# NOTE: STI uses 125 Hz .. 8 kHz, so the 63 Hz column is carried for
# completeness but is not used by the STI calculation.
NC_BAND_FREQUENCIES_HZ = [63, 125, 250, 500, 1000, 2000, 4000, 8000]

NC_CURVES = {
    15: [47, 36, 29, 22, 17, 14, 12, 11],
    20: [51, 40, 33, 26, 22, 19, 17, 16],
    25: [54, 44, 37, 31, 27, 24, 22, 21],
    30: [57, 48, 41, 35, 31, 29, 28, 27],
    35: [60, 52, 45, 40, 36, 34, 33, 32],
    40: [64, 56, 50, 45, 41, 39, 38, 37],
    45: [67, 60, 54, 49, 46, 44, 43, 42],
    50: [71, 64, 58, 54, 51, 49, 48, 47],
    55: [74, 67, 62, 58, 56, 54, 53, 52],
    60: [77, 71, 67, 63, 61, 59, 58, 57],
    65: [80, 75, 71, 68, 66, 64, 63, 62],
    70: [83, 79, 75, 72, 71, 70, 69, 68],
}


def nc_curve_levels(nc_value):
    """Returns {octave_hz: SPL_dB} for the seven STI octave bands
    (125 Hz - 8 kHz) of the given NC curve."""
    if nc_value not in NC_CURVES:
        raise ValueError("No NC curve defined for NC-%s" % nc_value)
    row = NC_CURVES[nc_value]
    by_band = dict(zip(NC_BAND_FREQUENCIES_HZ, row))
    return {hz: float(by_band[hz]) for hz in OCTAVE_BANDS_HZ}


def available_nc_values():
    return sorted(NC_CURVES.keys())


# ---------------------------------------------------------------------------
# Empirical %ALcons <-> STI conversion.
#
# PROVENANCE, stated plainly because it is weaker than the rest of this file:
# this is an empirical fit attributed to Farrel Becker. It is NOT part of
# IEC 60268-16 and no primary paper could be located for it. The publisher
# that carries it adds the caveat "this conversion is to be regarded with
# some doubts, because there are different evaluations of measurements".
#
# What CAN be shown is that the curve reproduces the published
# STI <-> %ALcons correspondence table to within ~0.6 percentage points:
#     STI 0.30 -> 33 % published, 33.56 % here
#     STI 0.45 -> 15 % published, 14.89 % here
#     STI 0.60 ->  7 % published,  6.60 % here
#     STI 0.75 ->  3 % published,  2.93 % here
#
# The two published forms are algebraically identical:
#     %ALcons = 170.5405 * exp(-5.419 * STI)
#     STI     = 0.9482 - 0.1845 * ln(%ALcons)
#
# Because of the weak provenance this is reported only as a SECONDARY
# cross-check. The primary %ALcons stays the energy-based Peutz / Bistafa &
# Bradley form, which is computed from the model's own physics.
BECKER_A = 170.5405
BECKER_B = 5.419

# Below about STI 0.2 the fit is extrapolation: it runs to 170 % where
# %ALcons is conventionally capped at 100 %.
BECKER_MIN_TRUSTED_STI = 0.20
ALCONS_MAX_PERCENT = 100.0


def alcons_from_sti(sti):
    """Empirical (Becker) %ALcons from STI. Returns (percent, trusted)."""
    value = BECKER_A * math.exp(-BECKER_B * float(sti))
    if value > ALCONS_MAX_PERCENT:
        value = ALCONS_MAX_PERCENT
    return value, (sti >= BECKER_MIN_TRUSTED_STI)


def sti_from_alcons(alcons_percent):
    """Inverse of alcons_from_sti; exact analytic round-trip."""
    if alcons_percent <= 0:
        return None
    return 0.9482 - 0.1845 * math.log(float(alcons_percent))


# ---------------------------------------------------------------------------
# Auditory masking and the absolute speech reception threshold.
#
# IEC 60268-16 models two hearing effects as imaginary masking noise that
# lowers the effective modulation index:
#
#   1. AUDITORY SPREAD OF MASKING -- a strong lower octave band masks the
#      band above it. Only the immediately lower band is considered, so the
#      lowest band (125 Hz) is never masked.
#   2. ABSOLUTE SPEECH RECEPTION THRESHOLD -- the hearing floor, a lower
#      limit on the effective noise in each band.
#
#         I_am,k = I_(k-1) * 10^(La/10)        La from the masker's LEVEL
#         I_rt,k = 10^(A_k/10)
#         m'     = m * I_k / (I_k + I_noise,k + I_am,k + I_rt,k)
#
# The noise term is the same one as before: I_k/(I_k+I_noise,k) is
# identically 1/(1+10^(-SNR/10)), so switching masking off reproduces
# apply_noise_snr() exactly.
#
# PROVENANCE. The piecewise-linear masking slope below is the standard's
# own specification -- IEC 60268-16 Annex A, Table A.2 -- not a fit derived
# by any implementer. The reception thresholds are Table A.3. Since the
# standard is paywalled, both were confirmed against accessible sources:
#
#   [1] H. Steeneken (co-originator of the STI), "Basics of the STI
#       measuring method", Table I (octave level specific slope of masking,
#       with amf = 10^(slope/10)) and Table II (absolute reception
#       threshold). This is the BANDED empirical data the standard's
#       continuous form smooths: -40, -35, -25, -20, -15, -10 dB/oct.
#   [2] Architectural Institute of Japan, benchmark problem A-21
#       (news-sv.aij.or.jp/kankyo/s24/benchmark/a21/a21_j.html), which
#       reproduces Table A.2 and Table A.3 verbatim, attributed to
#       IEC 60268-16:2020, alongside a full worked example.
#
# The two agree, and the piecewise form reproduces [1]'s banded table at
# every band centre (50 dB -> -40, 60 -> -35, 70 -> -24.8, 80 -> -19.8,
# 90 -> -14.8, >=100 -> -10). Note that [1] steps -35 straight to -25,
# skipping -30; the narrow 1.8*L - 146.9 segment, active only between 63
# and 67 dB, is what bridges that gap while staying continuous at both
# joins. That reconciliation is the standard's, which is why implementing
# it here carries no third-party licensing question.
#
# VALIDATED end to end against [2]: this implementation reproduces its
# published STI of 0.311 as 0.3107. See tests/BENCHMARK_A21.md.
ABSOLUTE_RECEPTION_THRESHOLD_DB = [46.0, 27.0, 12.0, 6.5, 7.5, 8.0, 12.0]


def masking_slope_db(level_db):
    """La: the octave-level-specific slope of masking, in dB per octave.

    Depends only on the LEVEL of the masking band, not on its frequency."""
    if level_db < 63.0:
        return 0.5 * level_db - 65.0
    if level_db < 67.0:
        return 1.8 * level_db - 146.9
    if level_db < 100.0:
        return 0.5 * level_db - 59.8
    return -10.0


def auditory_masking_factor(level_db):
    """amf = 10^(La/10): the intensity fraction a band of this level
    imposes on the octave band above it."""
    return 10.0 ** (masking_slope_db(level_db) / 10.0)


def modulation_reduction_factors(speech_spl_by_band, noise_spl_by_band=None,
                                 use_masking=True):
    """{hz: factor} such that m'(F) = m(F) * factor, per octave band.

    speech_spl_by_band must be ABSOLUTE band levels in dB SPL at the
    listener -- masking and the reception threshold are level dependent,
    so a relative scale would silently give the wrong answer.
    """
    noise_spl_by_band = noise_spl_by_band or {}
    factors = {}
    for index, hz in enumerate(OCTAVE_BANDS_HZ):
        speech_db = speech_spl_by_band.get(hz)
        if speech_db is None:
            factors[hz] = 1.0
            continue
        i_signal = 10.0 ** (speech_db / 10.0)

        i_noise = 0.0
        if hz in noise_spl_by_band:
            i_noise = 10.0 ** (noise_spl_by_band[hz] / 10.0)

        i_masking = 0.0
        i_threshold = 0.0
        if use_masking:
            if index > 0:                      # 125 Hz has no band below it
                lower_db = speech_spl_by_band.get(OCTAVE_BANDS_HZ[index - 1])
                if lower_db is not None:
                    i_masking = (10.0 ** (lower_db / 10.0)
                                 * auditory_masking_factor(lower_db))
            i_threshold = 10.0 ** (
                ABSOLUTE_RECEPTION_THRESHOLD_DB[index] / 10.0)

        total = i_signal + i_noise + i_masking + i_threshold
        factors[hz] = (i_signal / total) if total > 0 else 1.0
    return factors
