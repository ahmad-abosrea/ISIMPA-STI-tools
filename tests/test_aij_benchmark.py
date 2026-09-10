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
"""Validation against a PUBLISHED, third-party STI benchmark.

Architectural Institute of Japan, benchmark problem A-21, computed to
IEC 60268-16:2020 and agreed by three independent contributors (Hoshi,
Okubo, Nishikawa):

    https://news-sv.aij.or.jp/kankyo/s24/benchmark/a21/a21_j.html

    octave band (Hz)      125    250    500     1k     2k     4k     8k
    level incl. noise      50     45     42     40     36     36     30
    background noise       48     40     34     30     27     25     23
    published SNR       -2.33   3.35   7.25   9.54   8.42  10.64   6.03
    published MTI       0.131  0.284  0.330  0.324  0.322  0.319  0.292
    published STI                                                  0.311

These tests cover everything that does NOT need the 1.1 MB impulse
response: the STI assembly, the SNR convention, and the masking /
reception-threshold wiring. The full end-to-end run, which additionally
exercises octave filtering and the MTF integral from the raw impulse
response, reproduced STI = 0.3107 against the published 0.311 (worst
per-band MTI error 0.010, at 125 Hz). See BENCHMARK_A21.md.
"""
import math

import pytest

import sti_math as sm

BANDS = sm.OCTAVE_BANDS_HZ
LEVEL_WITH_NOISE = [50.0, 45.0, 42.0, 40.0, 36.0, 36.0, 30.0]
NOISE = [48.0, 40.0, 34.0, 30.0, 27.0, 25.0, 23.0]
PUBLISHED_SNR = [-2.33, 3.35, 7.25, 9.54, 8.42, 10.64, 6.03]
PUBLISHED_MTI = [0.131, 0.284, 0.330, 0.324, 0.322, 0.319, 0.292]
PUBLISHED_STI = 0.311


def _signal_only(level_with_noise, noise):
    """The published band levels INCLUDE the background noise."""
    return 10.0 * math.log10(10 ** (level_with_noise / 10.0)
                             - 10 ** (noise / 10.0))


def test_published_mti_reproduces_the_published_sti():
    """Validates the alpha/beta male weighting AND the square-root
    redundancy term against an external reference."""
    mti = dict(zip(BANDS, PUBLISHED_MTI))
    assert sm.compute_sti(mti) == pytest.approx(PUBLISHED_STI, abs=0.0005)


def test_the_redundancy_term_needs_its_square_root():
    """Without the square root the benchmark comes out at 0.388 instead of
    0.311 -- an error of 0.077, or 2.6 JND. Two secondary sources render
    the formula without the radical (it is drawn as vector graphics, not
    text, so extraction drops it); this test is why we do not."""
    m = list(PUBLISHED_MTI)
    no_sqrt = (sum(a * x for a, x in zip(sm.ALPHA_MALE, m))
               - sum(b * m[i] * m[i + 1] for i, b in enumerate(sm.BETA_MALE)))
    assert no_sqrt == pytest.approx(0.3878, abs=0.001)
    assert abs(no_sqrt - PUBLISHED_STI) > 20 * sm.STI_JND / 10.0


@pytest.mark.parametrize("index", range(7))
def test_snr_convention_matches_the_published_values(index):
    """Confirms both snr_from_levels() and the reading that the published
    levels include the noise."""
    signal = _signal_only(LEVEL_WITH_NOISE[index], NOISE[index])
    got = sm.snr_from_levels(signal, NOISE[index])
    assert got == pytest.approx(PUBLISHED_SNR[index], abs=0.01)


def test_reception_threshold_dominates_the_125_hz_band():
    """The benchmark's 125 Hz speech level (45.67 dB) sits essentially ON
    the 46 dB reception threshold, which is why that band is the most
    sensitive one in the whole calculation -- and why it carries the
    largest residual in the end-to-end run."""
    speech = {hz: _signal_only(l, n)
              for hz, l, n in zip(BANDS, LEVEL_WITH_NOISE, NOISE)}
    assert speech[125] == pytest.approx(45.67, abs=0.01)
    assert sm.ABSOLUTE_RECEPTION_THRESHOLD_DB[0] == 46.0

    i_sig = 10 ** (speech[125] / 10.0)
    i_thr = 10 ** (sm.ABSOLUTE_RECEPTION_THRESHOLD_DB[0] / 10.0)
    assert i_thr > i_sig                      # the floor outweighs the signal


def test_masking_and_threshold_change_the_answer_materially():
    """Guards against the hearing model being silently bypassed."""
    speech = {hz: _signal_only(l, n)
              for hz, l, n in zip(BANDS, LEVEL_WITH_NOISE, NOISE)}
    noise = dict(zip(BANDS, NOISE))

    with_model = sm.modulation_reduction_factors(speech, noise, use_masking=True)
    without = sm.modulation_reduction_factors(speech, noise, use_masking=False)

    assert with_model[125] < without[125]
    assert with_model[125] == pytest.approx(0.264, abs=0.002)
    for hz in BANDS:
        assert 0.0 < with_model[hz] <= without[hz] <= 1.0
