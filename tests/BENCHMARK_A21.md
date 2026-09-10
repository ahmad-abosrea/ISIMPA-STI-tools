# Validation against AIJ benchmark A-21

Third-party validation of this STI implementation against a **published
reference**, computed to IEC 60268-16:2020 and agreed by three independent
contributors (Hoshi, Okubo, Nishikawa).

Source: Architectural Institute of Japan, benchmark problem A-21
<https://news-sv.aij.or.jp/kankyo/s24/benchmark/a21/a21_j.html>

## Result

Run `python benchmark_a21_full.py` (needs `a21.fb` from the page above,
288000 float32 samples at 48 kHz, ~1.1 MB, placed alongside the script).

| octave band | published MTI | this implementation | delta |
|---|---|---|---|
| 125 Hz | 0.131 | 0.141 | +0.010 |
| 250 Hz | 0.284 | 0.283 | −0.001 |
| 500 Hz | 0.330 | 0.327 | −0.003 |
| 1 kHz | 0.324 | 0.323 | −0.001 |
| 2 kHz | 0.322 | 0.324 | +0.002 |
| 4 kHz | 0.319 | 0.318 | −0.001 |
| 8 kHz | 0.292 | 0.292 | +0.000 |
| **STI** | **0.311** | **0.3107** | **−0.0003** |

The STI agrees to within 0.0003 — about **1 % of one JND** (0.03).

## What this validates

The run drives the production function `sti_tool.ComputeSTI_SPPS()`, so
every stage of the real chain is exercised end to end:

- octave-band filtering of the impulse response
- the MTF integral at all 14 modulation frequencies
- background noise folding
- **auditory masking**
- **absolute speech reception threshold**
- male alpha/beta weighting with the square-root redundancy term
- MTI averaging and final STI assembly

## The 125 Hz residual

The largest per-band deviation, +0.010, is at 125 Hz, and that band is
expected to be the least stable one in this particular benchmark:

| band | I_signal | I_noise | I_threshold | reduction factor |
|---|---|---|---|---|
| 125 Hz | 26 % | 45 % | **28 %** | 0.264 |
| 250 Hz | 67 % | 31 % | 2 % | 0.673 |
| 1 kHz | 90 % | 10 % | 0 % | 0.900 |

The benchmark's 125 Hz speech level works out at **45.67 dB**, essentially
identical to the **46 dB** absolute reception threshold. The hearing floor
therefore contributes more intensity than the speech does, and the band
sits on a knife edge: a small difference in band energy moves the result
noticeably. Every other band has the threshold contributing 0–2 %.

The most likely cause of that small energy difference is the octave filter.
This implementation isolates bands by zeroing out-of-band FFT bins (edges
at `fc/sqrt(2)` .. `fc*sqrt(2)`), a brick-wall approximation, whereas the
standard assumes an IEC 61260 filter with finite skirts that lets a little
adjacent-band energy through. At 125 Hz, where the level sits on the
threshold, that difference is amplified; elsewhere it is negligible, which
is exactly the pattern observed.

## Permanent regression tests

`test_aij_benchmark.py` pins the parts that need no audio file:

- the published MTI must reproduce the published STI (validates alpha,
  beta, and the square root)
- dropping the square root must give 0.388, not 0.311 — an error of 0.077,
  or 2.6 JND
- the published per-band SNR must be reproduced, which also confirms that
  the published band levels *include* the background noise
- the hearing model must materially change the answer, so it cannot be
  silently bypassed
