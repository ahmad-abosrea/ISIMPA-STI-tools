# ISIMPA-STI-tools

[![tests](https://github.com/ahmad-abosrea/ISIMPA-STI-tools/actions/workflows/tests.yml/badge.svg)](https://github.com/ahmad-abosrea/ISIMPA-STI-tools/actions/workflows/tests.yml)

Speech-intelligibility post-processing for I-Simpa.

Python plugins that add **STI** (Speech Transmission Index, IEC 60268-16)
and **%ALcons** (articulation loss of consonants) to
[I-Simpa](https://i-simpa.univ-gustave-eiffel.fr), for punctual receivers
and as a colour map on surface receivers.

Validated against a published third-party benchmark: it reproduces the
Architectural Institute of Japan's benchmark A-21 STI of **0.311** as
**0.3107** — about 1 % of one JND. See [tests/BENCHMARK_A21.md](tests/BENCHMARK_A21.md).

---

## What you get

Right-click a receiver in the **Results** tree:

| Tool | Menu entry | Where |
|---|---|---|
| `sti_tool` | **Compute STI** | each punctual receiver |
| `sti_tool` | **Compute STI Map** | a surface-receiver folder |
| `alcons_tool` | **Compute %ALcons** | each punctual receiver |
| `debug_tool` | **DEBUG: Dump …** | receivers and scene sources |

`debug_tool` is a diagnostic aid: it prints the real structure of I-Simpa's
result tree and data files. It is not needed for normal use, but it is what
makes it possible to check assumptions against the live application rather
than guessing.

## Scope — read this before quoting a number

- **Full STI (STI-14)**: 7 octave bands × 14 modulation frequencies.
  **Not** STIPA, RASTI or STITEL.
- **Male weighting only.** IEC specifies male coefficients for this purpose.
- **Indirect method**: the modulation transfer function is obtained from the
  impulse response, not by playing a test signal through the system.
- **Includes** background noise, auditory masking, and the absolute speech
  reception threshold.
- **Two engines**: for SPPS results the real echogram is integrated; for TCR
  results the analytical exponential MTF is used from the per-band T60. For
  SPPS a second estimate from the Schroeder decay is reported alongside, and
  the difference between them is diagnostic.

This is an independent implementation. It is **not** endorsed, certified or
verified by IEC. Anyone needing normative certainty must consult
IEC 60268-16 itself.

## Install

Copy the tool folders into I-Simpa's `UserScript` directory:

```
UserScript/sti_tool      ->  <I-Simpa>/UserScript/sti_tool
UserScript/alcons_tool   ->  <I-Simpa>/UserScript/alcons_tool
UserScript/debug_tool    ->  <I-Simpa>/UserScript/debug_tool
```

On Windows that is typically `C:\Program Files\I-SIMPA\UserScript\`, which
usually needs administrator rights. **Restart I-Simpa** afterwards — the
startup scripts are only loaded at launch.

Copy whole folders. `sti_tool` in particular needs all five files;
`alcons_tool` imports `sti_noise` and `sti_math` from it for the STI
cross-check, and degrades gracefully if `sti_tool` is absent.

## Getting trustworthy numbers out of it

Two model settings dominate the result, and both bite quietly:

**Time step.** %ALcons splits the echogram into direct and reverberant
energy, and one time step per source is counted as direct. With `N` sources
the direct window is `N × dt`, and that must stay small against the decay
constant `tau = T60/13.8`. A workable target is `N × dt <= 0.2·tau`. With 11
sources, a 0.7 s T60 and a 10 ms step, 110 ms gets labelled "direct" —
around 90 % of all the energy — and %ALcons comes out far too optimistic.
At 1 ms the same case is usable.

**Simulation length.** With a particle-extinction limit of 10^n, the
echogram cannot carry more than `n × 10` dB of decay however long you run.
Time spent past that is empty bins. Spend the budget on resolution instead.

STI is far less sensitive to both, because the MTF integral never needs a
direct/reverberant split. **When the two %ALcons estimates diverge by more
than a factor of two, the tool says so**, and in a multi-source model the
STI-derived value is the one to trust.

## Tests

```
python -m pytest tests/ -q
```

145 tests, no third-party dependencies. They are run automatically on
Python 3.8 (the version I-Simpa embeds) and on a current Python. They run against fake `uictrl` and
`libsimpa` modules that reproduce I-Simpa's real data layout, so no
installation of I-Simpa is needed.

`tests/benchmark_a21_full.py` runs the full published benchmark end to end.
It needs `a21.fb` from the AIJ page linked in `tests/BENCHMARK_A21.md`; that
file is not redistributed here.

## Licence

**GPL-3.0-or-later** — see [LICENSE](LICENSE).

These tools run inside I-Simpa and use its API, and I-Simpa is
GPL-3.0-or-later, so this licence is required rather than chosen.

Standards constants, formulas, attributions and the reasoning behind each
are recorded in [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md). No
third-party source code is included in this project.
