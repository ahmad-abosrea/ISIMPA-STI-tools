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
import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _userscript_dir():
    return os.path.join(_REPO, "UserScript")


def _tool_dir(name):
    return os.path.join(_REPO, "UserScript", name)


"""AIJ benchmark A-21 -- full end-to-end validation of the STI chain.

Published reference (Architectural Institute of Japan, conforming to
IEC 60268-16:2020), agreed by three independent contributors:

    MTI  125..8k : 0.131 0.284 0.330 0.324 0.322 0.319 0.292
    STI          : 0.311

Input: a21.fb, 288000 float32 samples at 48 kHz, mono.
Band levels INCLUDE the background noise, so the signal-only level is
recovered as 10*log10(10^(L/10) - 10^(N/10)) -- a convention already
confirmed by reproducing their published per-band SNR exactly.

Octave filtering is done by zeroing out-of-band FFT bins (edges at
fc/sqrt(2) .. fc*sqrt(2)). That is a brick-wall approximation to the
IEC 61260 filter the standard assumes; if a band disagrees while the
others match, the filter is the first thing to suspect.
"""
import array
import cmath
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import isimpa_fixtures as fx
fx.install()
sys.path.insert(0, _tool_dir('sti_tool'))
sys.path.insert(0, _userscript_dir())
import sti_math as sm
import sti_tool                                    # the production module

FS = 48000.0
PUBLISHED_MTI = [0.131, 0.284, 0.330, 0.324, 0.322, 0.319, 0.292]
PUBLISHED_STI = 0.311
LEVEL_WITH_NOISE = [50.0, 45.0, 42.0, 40.0, 36.0, 36.0, 30.0]
NOISE = [48.0, 40.0, 34.0, 30.0, 27.0, 25.0, 23.0]
ENVELOPE_BIN_S = 0.001          # 1 ms: far finer than the 12.5 Hz max


def load_ir(path):
    raw = open(path, "rb").read()
    a = array.array("f")
    a.frombytes(raw)
    if sys.byteorder != "little":
        a.byteswap()
    return list(a)


def fft(x, inverse=False):
    """Iterative in-place radix-2 FFT. len(x) must be a power of two."""
    n = len(x)
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j |= bit
        if i < j:
            x[i], x[j] = x[j], x[i]
    length = 2
    while length <= n:
        ang = (2 if inverse else -2) * math.pi / length
        wl = cmath.exp(complex(0, ang))
        half = length >> 1
        for i in range(0, n, length):
            w = complex(1, 0)
            for k in range(i, i + half):
                u = x[k]
                v = x[k + half] * w
                x[k] = u + v
                x[k + half] = u - v
                w *= wl
        length <<= 1
    if inverse:
        for i in range(n):
            x[i] /= n
    return x


def main():
    t0 = time.time()
    ir = load_ir(os.path.join(HERE, "a21.fb"))
    n_src = len(ir)
    size = 1
    while size < n_src:
        size <<= 1
    print("IR: %d samples, FFT size %d" % (n_src, size))

    spec = [complex(v, 0.0) for v in ir] + [0j] * (size - n_src)
    fft(spec)
    print("forward FFT done (%.1fs)" % (time.time() - t0))

    df = FS / size
    echogram = {}
    for idx, hz in enumerate(sm.OCTAVE_BANDS_HZ):
        lo, hi = hz / math.sqrt(2.0), hz * math.sqrt(2.0)
        k_lo, k_hi = int(math.floor(lo / df)), int(math.ceil(hi / df))
        band = [0j] * size
        for k in range(k_lo, min(k_hi + 1, size // 2)):
            band[k] = spec[k]
            if k:
                band[size - k] = spec[size - k]        # keep it real
        fft(band, inverse=True)

        # squared band-limited IR, averaged into 1 ms energy bins
        step = int(round(ENVELOPE_BIN_S * FS))
        times, h2 = [], []
        for start in range(0, n_src, step):
            acc = 0.0
            for i in range(start, min(start + step, n_src)):
                r = band[i].real
                acc += r * r
            times.append((start + step * 0.5) / FS)
            h2.append(acc)
        echogram[hz] = (times, h2)
        print("  %5d Hz band done (%.1fs)" % (hz, time.time() - t0))

    speech = {}
    noise = {}
    for i, hz in enumerate(sm.OCTAVE_BANDS_HZ):
        sig = 10.0 * math.log10(10 ** (LEVEL_WITH_NOISE[i] / 10.0)
                                - 10 ** (NOISE[i] / 10.0))
        speech[hz] = sig
        noise[hz] = NOISE[i]

    factors = sm.modulation_reduction_factors(speech, noise, use_masking=True)
    sti, mti = sti_tool.ComputeSTI_SPPS(echogram, factors)   # production code

    print()
    print("=== AIJ benchmark A-21 ===")
    print("  band      published      ours      delta")
    worst = 0.0
    for i, hz in enumerate(sm.OCTAVE_BANDS_HZ):
        d = mti[hz] - PUBLISHED_MTI[i]
        worst = max(worst, abs(d))
        print("  %5d Hz     %.3f      %.3f    %+.3f" % (hz, PUBLISHED_MTI[i], mti[hz], d))
    print()
    print("  STI          %.3f      %.4f    %+.4f" % (PUBLISHED_STI, sti, sti - PUBLISHED_STI))
    print("  worst per-band MTI error: %.3f" % worst)
    print("  JND for STI is 0.03; elapsed %.1fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
