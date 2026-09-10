# Third-party notices and attributions

This project contains **no third-party source code**. Everything here was
written for this project. What follows records the published standards,
formulas and reference data the implementation relies on, and why each is
used the way it is.

---

## 1. I-Simpa — the host application (GPL-3.0-or-later)

These tools are plugins for **I-Simpa**, a 3D sound-propagation modelling
GUI by UMRAE, CEREMA and Université Gustave Eiffel
(<https://i-simpa.univ-gustave-eiffel.fr>), licensed **GPL-3.0-or-later**.

They run inside I-Simpa's embedded Python interpreter and call its
`uictrl` and `libsimpa` APIs, and the menu-registration boilerplate
(`class manager` / `getmenu` / `register_menu_manager`) follows the plugin
interface demonstrated by I-Simpa's own bundled example scripts.

**This project is therefore licensed GPL-3.0-or-later**, and its `LICENSE`
file is the same GNU GPL v3 text that I-Simpa ships. This is a requirement
of that licence, not a preference.

I-Simpa is not affiliated with this project and does not endorse it.

---

## 2. IEC 60268-16 — Speech Transmission Index

The STI implementation follows **IEC 60268-16** ("Sound system equipment –
Part 16: Objective rating of speech intelligibility by speech transmission
index"), Edition 4:2011 / Edition 5:2020.

The following **numeric constants** from the standard are reproduced,
because interoperating with a standard is impossible without them:

| Quantity | Source in the standard |
|---|---|
| Male octave weighting α (7 values) and redundancy β (6 values) | Table A.1 |
| Slope of masking, piecewise-linear in the masker level | Table A.2 |
| Absolute speech reception threshold, 7 values | Table A.3 |
| The seven octave bands and fourteen modulation frequencies | normative text |

**No text, tables as laid out, figures or other expressive content from the
standard is reproduced.** Only the bare numeric values needed to compute
the index are used, together with the formulas, which are mathematical
facts rather than expression.

IEC 60268-16 is copyright IEC and must be purchased from IEC or a national
standards body. **Anyone needing normative certainty must consult the
standard itself** — this project is an independent implementation and is
not endorsed, certified or verified by IEC.

Scope actually implemented: **full STI (STI-14), male weighting, indirect
method** (MTF derived from an impulse response). Not STIPA, not RASTI, not
STITEL, and no female weighting.

---

## 3. Corroborating sources for the IEC constants

Because the standard is paywalled, every constant was checked against at
least two independent published sources before use:

- **H. J. M. Steeneken**, *"Basics of the STI measuring method"*
  (<https://www.steeneken.nl>) — Table I (octave-level-specific slope of
  masking, and the corresponding masking factors), Table II (reception
  threshold, and the male/female weighting factors). Steeneken co-originated
  the STI method.
- **Architectural Institute of Japan**, benchmark problem **A-21**
  (<https://news-sv.aij.or.jp/kankyo/s24/benchmark/a21/a21_j.html>) —
  reproduces Tables A.1–A.4 and a complete worked example conforming to
  IEC 60268-16:2020.
- Published Edition 5 compliance documentation accompanying a widely used
  numerical-computing environment — α, β, the octave bands and the
  modulation frequency generator `10^((-2:11)/10)`.

### A note on `rajmic/Speech-Transmission-Index-STI`

During development, the piecewise-linear masking function was
**cross-checked against** the reference implementation at
<https://github.com/rajmic/Speech-Transmission-Index-STI>
(© Šimon Cieslar, Pavel Záviška, Jiří Schimmel, Pavel Rajmic, Brno
University of Technology, 2023–2025).

That repository carries **no licence file**, so no code from it may be
reused, and **none is used here**. It is recorded only because honest
provenance requires it. The formula itself is not theirs to license: their
own code cites it as IEC Annex A, Table A.2, and the Architectural
Institute of Japan benchmark independently reproduces the same piecewise
function attributed to the same table. It is the standard's specification,
and this implementation follows the standard.

---

## 4. ISO 9921 — intelligibility rating scale

The five-step rating (bad / poor / fair / good / excellent) with boundaries
at STI 0.30, 0.45, 0.60 and 0.75, and the associated phonetically balanced
word-score ranges, follow **ISO 9921** ("Ergonomics – Assessment of speech
communication"). Numeric boundaries only.

The lettered "qualification band" scale (A+ … U) that appeared in earlier
development versions was **removed**: its thresholds could not be verified
against the normative text, and presenting it as IEC nomenclature claimed
more authority than could be supported.

---

## 5. NC curves

Octave-band noise-criterion limits (NC 15 … NC 70) follow the curves
originally due to **L. L. Beranek**, as tabulated in **ANSI/ASA S12.2** and
reproduced in standard HVAC-acoustics references. Numeric values only.

---

## 6. %ALcons — articulation loss of consonants

- The **energetic form**, `%ALcons = 8.9·T60·E_rev / (13.82·E_dir)`, follows
  Bies & Hansen, *Engineering Noise Control*, eq. (7.130), after
  **Bistafa & Bradley (2000)**.
- The **direct-to-reverberant rating scale**, eq. (7.131), is from the same
  source.
- The **9·T60 saturation ceiling** is **Peutz's** "without direct sound"
  limit.
- The **STI → %ALcons conversion**, `%ALcons = 170.5405·e^(−5.419·STI)`
  (equivalently `STI = 0.9482 − 0.1845·ln %ALcons`), is an **empirical
  relation attributed to Farrel Becker**. It is **not** part of
  IEC 60268-16, no primary publication for it could be located, and the
  source that carries it notes that it "is to be regarded with some doubts,
  because there are different evaluations of measurements". It is therefore
  reported only as a clearly labelled secondary cross-check.

Formulas are mathematical facts and are used with attribution; no text from
any of these works is reproduced.

---

## 7. Validation reference

The implementation is validated end to end against **Architectural
Institute of Japan benchmark problem A-21**, a published test case computed
to IEC 60268-16:2020 and agreed by three independent contributors. This
implementation reproduces the published STI of 0.311 as **0.3107**. See
`tests/BENCHMARK_A21.md`.

The benchmark's impulse-response files are **not redistributed** with this
project; the test that needs them downloads nothing and is skipped unless
you place the file yourself.

---

## 8. Trademarks

I-Simpa, IEC, ISO, ANSI/ASA and any product names mentioned are the
property of their respective owners, and are used here descriptively, to
identify the software or standard being referred to. No affiliation or
endorsement is implied.
