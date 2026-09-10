# -*- coding: cp1252 -*-
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
# STI (Speech Transmission Index, IEC 60268-16 full method) post-processing
# tool. Adds a "Compute STI" entry to the right-click menu of each
# INDIVIDUAL punctual receiver folder (e.g. "Receiver 1").
#
# Branches on where the results came from:
#   - SPPS (particle tracing): uses the REAL echogram (h^2(t) per octave
#     band, read from the receiver's "Sound level" report) and numerically
#     integrates the true MTF.
#   - TCR (classical reverberation theory, no real impulse response): uses
#     the analytical exponential-decay MTF from the per-band T60.
#
# GRID LAYOUT (verified live against I-Simpa 2.x via a debug dump, NOT
# assumed -- an earlier version of this file had the axes the wrong way
# round and silently found no data):
#
#   'Sound level' report:
#     row[0] = ['', '10.0 ms', '20.0 ms', ..., '2000.0 ms', 'Total']
#     row[i] = ['63 Hz', <one value per time step...>, <total>]
#     => ROWS are frequency bands, COLUMNS are time steps.
#     Values are dB; converted back to linear energy (10**(dB/10)) here
#     because the MTF integral needs a quantity homogeneous to p^2.
#     The trailing 'Total' column and the 'Global'/'Average' rows are
#     cumulative summaries and are excluded.
#
#   'Acoustic parameters' report:
#     row[0] = ['', 'Sound level (dB)', ..., 'RT-30 (s)', 'EDT (s)', ...]
#     row[i] = ['63 Hz', <values...>]
#     Note the RT column is 'RT-30 (s)' in this build while the shipped
#     sample script uses the translated 'TR-%g (s)'; both are accepted.
#
# The tool triggers I-Simpa's own acoustic-parameter calculation itself for
# the TCR branch, so RT does not have to be computed by hand first.
#
# The heavy math (IEC 60268-16 alpha/beta weighting, both MTF formulas,
# noise folding, STI combination) lives in sti_math.py (a plain-Python,
# pytest-validated module alongside this file) -- see
# the tests/ folder of this repository.
#
# NOTE: ambient/electrical noise per band is not exposed by the uictrl
# Python API, so this tool asks for a flat SNR via a dialog.
from __future__ import print_function
import math
import os
import re
import sys
import uictrl as ui
from libsimpa import *

# sti_math.py lives alongside this file (inside the sti_tool/ package
# folder), but only the top-level UserScript/ folder is guaranteed to be on
# sys.path when I-Simpa's loader runs __ui_startup__.py -- a plain
# 'import sti_math' fails with ModuleNotFoundError otherwise (confirmed
# live). Adding this file's own directory to sys.path makes it robust.
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)
import sti_math as sm

SPPS_LABELS = (u"SPPS",)
TCR_LABELS = (u"TCR", u"Diffusion model")
PUNCTUAL_RECEIVERS_FOLDER_NAME = u"Punctual receivers"
# How far the nearest available band may be from an IEC target band before
# we refuse to use it (e.g. a project with no 8 kHz band should error, not
# silently reuse the 4 kHz result).
BAND_MATCH_TOLERANCE_RATIO = 0.30

MASKING_LABEL = u"Auditory masking + reception threshold"
MASKING_ON = u"On - full IEC hearing model (recommended)"
MASKING_OFF = u"Off - reverberation and ambient noise only"
MASKING_MODES = [MASKING_ON, MASKING_OFF]

NOISE_MODE_NONE = u"None - noise-free STI (reverberation only)"
NOISE_MODE_USER = u"User Defined"


def _nc_mode_label(nc):
    return u"NC %d" % nc


# Dropdown order: noise-free, user defined, then NC 15 .. NC 65.
# "None" is kept as an option because the reverberation-only result is a
# useful upper bound, and because it is the only way to see the room's
# own contribution separated from the noise.
NOISE_MODES = ([NOISE_MODE_NONE, NOISE_MODE_USER]
               + [_nc_mode_label(nc) for nc in sm.available_nc_values()])

DEFAULT_NOISE_SPL_DB = 30.0
# "noise-free" is expressed as a very large SNR rather than a special case,
# so the same code path runs either way: 10^(-99/10) is ~1e-10, i.e. the
# noise term in m'(F) = m(F)/(1+10^(-SNR/10)) vanishes.
NOISE_FREE_SNR_DB = 99.0


# The dialog sorts its fields alphabetically, so a plain "125 Hz" label
# would list 1000 before 125. Zero-padding the frequency makes the
# alphabetical order match the numeric order: 0063, 0125, 0250, 0500,
# 1000, 2000, 4000, 8000.
#
# 63 Hz is shown even though STI uses 125 Hz - 8 kHz only, so the entered
# spectrum matches the familiar 8-band NC table; the 63 Hz value is read
# and echoed back but never reaches the STI calculation.
NOISE_ENTRY_BANDS_HZ = [63] + list(sm.OCTAVE_BANDS_HZ)


def _noise_label(hz):
    return u"Noise SPL %04d Hz (dB)" % hz


_BAND_RE = re.compile(r"(\d+(?:[.,]\d+)?)\s*(k?)Hz", re.IGNORECASE)


def _band_label_to_hz(label):
    if not isinstance(label, (str, bytes)):
        return None
    if isinstance(label, bytes):
        try:
            label = label.decode("cp1252")
        except Exception:
            return None
    m = _BAND_RE.search(label)
    if not m:
        return None
    value = float(m.group(1).replace(",", "."))
    if m.group(2).lower() == "k":
        value *= 1000.0
    return value


def _time_label_to_seconds(label):
    """'10.0 ms' -> 0.010. Returns None for 'Total' and the empty corner
    cell, which is how those non-time columns get excluded."""
    if not isinstance(label, (str, bytes)):
        return None
    if isinstance(label, bytes):
        try:
            label = label.decode("cp1252")
        except Exception:
            return None
    s = label.strip()
    if not s.lower().endswith("ms"):
        return None
    try:
        return float(s[:-2].strip().replace(",", ".")) / 1000.0
    except ValueError:
        return None


def _db_to_linear_energy(value):
    """Converts an echogram cell (dB) to linear energy.

    -inf dB means literally zero energy in that step, which I-Simpa emits
    for steps before the direct sound arrives (confirmed in a live dump).
    That is VALID data meaning zero, so it maps to 0.0 and the time step
    is kept -- dropping those rows would leave a ragged time axis and
    distort the MTF integral's time weighting. NaN is a real failure."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f != f:
        return None
    if f == float("-inf"):
        return 0.0
    if f == float("inf"):
        return None
    return 10 ** (f / 10.0)


def receiver_output_path(receiver_folder_id, filename):
    """buildfullpath() lives on e_file, NOT on the plain element wrapper
    -- calling it on element raises AttributeError (confirmed live)."""
    return ui.e_file(receiver_folder_id).buildfullpath() + filename


def _as_finite_float(value):
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f != f:  # NaN
        return None
    if f == float("inf") or f == float("-inf"):
        return None
    return f


def _find_rt_column(header, tr_value=30):
    candidates = []
    try:
        candidates.append(ui._("TR-%g (s)") % tr_value)
    except Exception:
        pass
    candidates.append("RT-%g (s)" % tr_value)
    candidates.append("TR-%g (s)" % tr_value)
    for cand in candidates:
        if cand in header:
            return header.index(cand)
    pattern = re.compile(r"(?:RT|TR)\s*-?\s*%d\b" % tr_value, re.IGNORECASE)
    for idx, cell in enumerate(header):
        if isinstance(cell, (str, bytes)) and pattern.search(
            cell if isinstance(cell, str) else cell.decode("cp1252", "replace")
        ):
            return idx
    return -1


def DetectCalculationMode(element_id):
    """Walks up the ancestor chain looking for an SPPS or TCR/Diffusion
    model ancestor label."""
    el = ui.element(element_id)
    while True:
        infos = el.getinfos()
        label = infos.get("label", infos.get("name", ""))
        if label in SPPS_LABELS:
            return "SPPS"
        if label in TCR_LABELS:
            return "TCR"
        parentid = infos.get("parentid", -1)
        if parentid is None or parentid < 0:
            return None
        el = ui.element(parentid)


def FindSoundLevelGabeId(receiver_folder_id):
    """The COMBINED echogram (all sources summed) for this receiver.

    Deliberately looks only at DIRECT children. "Echogram per source"
    makes SPPS write an extra report, also named "Sound level", inside a
    sub-folder per source (reportmanager.cpp:620), and I-Simpa's
    getallelementbytype() recurses (element.cpp:314) -- so a first-match
    search could return ONE source's echogram instead of the sum. That
    would be wrong in a way that is hard to see, because the direct-arrival
    logic would then mark every source's arrival in an echogram containing
    only one source's energy.
    """
    folder = ui.element(receiver_folder_id)
    try:
        children = folder.childs()
    except Exception:
        children = []
    for child in children:
        try:
            infos = ui.element(child[0]).getinfos()
        except Exception:
            continue
        if (infos.get("typeElement") == ui.element_type.ELEMENT_TYPE_REPORT_GABE_RECP
                and infos.get("name") == "Sound level"):
            return child[0]
    return -1


def CountPerSourceEchograms(receiver_folder_id):
    """How many extra per-source 'Sound level' reports sit BELOW this
    receiver. Non-zero means "Echogram per source" was ticked; the
    combined report is still used, and this is reported so the choice is
    visible rather than silent."""
    folder = ui.element(receiver_folder_id)
    try:
        every = folder.getallelementbytype(
            ui.element_type.ELEMENT_TYPE_REPORT_GABE_RECP)
    except Exception:
        return 0
    total = 0
    for gid in every:
        try:
            if ui.element(gid).getinfos().get("name") == "Sound level":
                total += 1
        except Exception:
            continue
    return max(0, total - 1)


def _match_band(available_hz, target_hz):
    """Nearest available band to an IEC target, or None if nothing is
    close enough (so a missing band is reported rather than faked)."""
    if not available_hz:
        return None
    best = min(available_hz, key=lambda hz: abs(hz - target_hz))
    if abs(best - target_hz) > BAND_MATCH_TOLERANCE_RATIO * target_hz:
        return None
    return best


def GetT60PerOctaveBand_TCR(receiver_folder_id, sound_level_gabe_id):
    """Returns (t60_by_iec_band, missing_bands, invalid_bands)."""
    ui.application.sendevent(
        ui.element(sound_level_gabe_id),
        ui.idevent.IDEVENT_RECP_COMPUTE_ACOUSTIC_PARAMETERS,
        {"TR": "30"},
    )
    receiver = ui.element(receiver_folder_id)
    params_id = receiver.getelementbylibelle("Acoustic parameters")
    if params_id == -1:
        return None, None, None
    grid = ui.application.getdataarray(ui.element(params_id))
    if not grid or len(grid) < 2:
        return None, None, None

    header = grid[0]
    idcol = _find_rt_column(header, 30)
    if idcol == -1:
        print(ui._("Could not find an RT-30 column. Header was: %r") % (header,))
        return None, None, None

    available = {}
    invalid_source_bands = set()
    for row in grid[1:]:
        hz = _band_label_to_hz(row[0])
        if hz is None or idcol >= len(row):
            continue
        t60 = _as_finite_float(row[idcol])
        if t60 is None or t60 <= 0:
            invalid_source_bands.add(hz)
        else:
            available[hz] = t60

    result = {}
    missing = []
    invalid = []
    for target in sm.OCTAVE_BANDS_HZ:
        matched = _match_band(list(available.keys()), target)
        if matched is None:
            # distinguish "band exists but NaN" from "band absent entirely"
            if _match_band(list(invalid_source_bands), target) is not None:
                invalid.append(target)
            else:
                missing.append(target)
            continue
        result[target] = available[matched]
    return result, missing, invalid


def GetEchogramPerOctaveBand_SPPS(sound_level_gabe_id):
    """Returns (echogram_by_iec_band, missing_bands, decay_db_by_band).
    echogram maps octave_hz -> (times_s, h_squared)."""
    grid = ui.application.getdataarray(ui.element(sound_level_gabe_id))
    if not grid or len(grid) < 2:
        return None, None, None

    header = grid[0]
    # COLUMNS are time steps (excluding the trailing 'Total' summary).
    time_cols = []
    for idx in range(1, len(header)):
        t = _time_label_to_seconds(header[idx])
        if t is not None:
            time_cols.append((t, idx))
    if len(time_cols) < 2:
        return None, None, None
    time_cols.sort(key=lambda pair: pair[0])

    # ROWS are frequency bands (excluding 'Global'/'Average' summaries).
    band_rows = {}
    for row in grid[1:]:
        hz = _band_label_to_hz(row[0])
        if hz is not None:
            band_rows[hz] = row

    result = {}
    missing = []
    decay_db = {}
    for target in sm.OCTAVE_BANDS_HZ:
        matched = _match_band(list(band_rows.keys()), target)
        if matched is None:
            missing.append(target)
            continue
        row = band_rows[matched]
        times = []
        h2 = []
        levels_db = []
        for t, idx in time_cols:
            if idx >= len(row):
                continue
            e = _db_to_linear_energy(row[idx])
            if e is None:
                continue
            times.append(t)
            h2.append(e)
            if e > 0:
                levels_db.append(10 * math.log10(e))
        if len(times) < 2 or sum(h2) <= 0 or not levels_db:
            missing.append(target)
            continue
        result[target] = (times, h2)
        decay_db[target] = max(levels_db) - levels_db[-1]
    return result, missing, decay_db


def GetReceiverSpeechLevelPerBand(sound_level_gabe_id):
    """Speech SPL at this receiver, per octave band, taken from the 'Total'
    column of the Sound level report (the energetic sum over all time steps
    = the steady-state level from all sources, with their directivities
    already applied).

    This is the 'signal' half of the SNR. Returns {octave_hz: spl_db}."""
    grid = ui.application.getdataarray(ui.element(sound_level_gabe_id))
    if not grid or len(grid) < 2:
        return None
    header = grid[0]

    total_col = None
    for idx in range(len(header) - 1, 0, -1):
        cell = header[idx]
        if isinstance(cell, bytes):
            try:
                cell = cell.decode("cp1252")
            except Exception:
                continue
        if isinstance(cell, str) and cell.strip().lower() in ("total", "cumul"):
            total_col = idx
            break
    if total_col is None:
        total_col = len(header) - 1  # the last column is the cumulative total

    available = {}
    for row in grid[1:]:
        hz = _band_label_to_hz(row[0])
        if hz is None or total_col >= len(row):
            continue
        v = _as_finite_float(row[total_col])
        if v is not None:
            available[hz] = v
    if not available:
        return None

    result = {}
    for target in sm.OCTAVE_BANDS_HZ:
        matched = _match_band(list(available.keys()), target)
        if matched is not None:
            result[target] = available[matched]
    return result


def ComputeSTI_TCR(t60_by_band, reduction_factors):
    """reduction_factors[hz] folds in ambient noise and, when enabled,
    auditory masking and the absolute reception threshold."""
    mti_by_band = {}
    for hz in sm.OCTAVE_BANDS_HZ:
        t60 = t60_by_band[hz]
        ti_vals = []
        for f_mod in sm.MODULATION_FREQS_HZ:
            m = sm.mtf_analytical_exponential_decay(t60, f_mod)
            m_prime = m * reduction_factors[hz]
            ti_vals.append(sm.transmission_index(m_prime))
        mti_by_band[hz] = sm.mti_for_band(ti_vals)
    return sm.compute_sti(mti_by_band), mti_by_band


def ComputeSTI_SPPS(echogram_by_band, reduction_factors):
    mti_by_band = {}
    for hz in sm.OCTAVE_BANDS_HZ:
        times, h2 = echogram_by_band[hz]
        ti_vals = []
        for f_mod in sm.MODULATION_FREQS_HZ:
            m = sm.mtf_from_echogram(times, h2, f_mod)
            m_prime = m * reduction_factors[hz]
            ti_vals.append(sm.transmission_index(m_prime))
        mti_by_band[hz] = sm.mti_for_band(ti_vals)
    return sm.compute_sti(mti_by_band), mti_by_band


def sti_for_receiver(receiver_folder_id, gid, mode, recv_name,
                     noise_spl_by_band, use_masking=True):
    """STI for one receiver, without any dialog or file output.

    Shared by ComputeSTI and by alcons_tool's empirical cross-check, so
    the two tools can never report a different STI for the same receiver.
    Returns a dict, or None if the data needed was not available (the
    reason is printed, exactly as it was when this lived inline).
    """
    use_noise = bool(noise_spl_by_band)

    speech_spl_by_band = None
    # Masking and the reception threshold are LEVEL dependent, so the
    # absolute band levels are needed even when there is no ambient noise.
    if use_noise or use_masking:
        speech_spl_by_band = GetReceiverSpeechLevelPerBand(gid)
        if speech_spl_by_band is None or len(speech_spl_by_band) < len(sm.OCTAVE_BANDS_HZ):
            if use_noise:
                print(ui._(
                    "Could not read the speech level per band at this receiver, so "
                    "the signal-to-noise ratio cannot be formed. STI was NOT "
                    "computed."))
                return None
            print(ui._(
                "   Could not read the speech level per band, so auditory masking "
                "and the reception threshold were SKIPPED (they are level "
                "dependent). The result is reverberation-only."))
            speech_spl_by_band = None
            use_masking = False

    if speech_spl_by_band:
        reduction_factors = sm.modulation_reduction_factors(
            speech_spl_by_band, noise_spl_by_band if use_noise else None,
            use_masking)
        snr_db_by_band = {
            hz: (sm.snr_from_levels(speech_spl_by_band[hz], noise_spl_by_band[hz])
                 if use_noise else NOISE_FREE_SNR_DB)
            for hz in sm.OCTAVE_BANDS_HZ
        }
    else:
        reduction_factors = {hz: 1.0 for hz in sm.OCTAVE_BANDS_HZ}
        snr_db_by_band = {hz: NOISE_FREE_SNR_DB for hz in sm.OCTAVE_BANDS_HZ}

    truncation_warning = None
    sti_schroeder = None
    t60_schroeder = None

    if mode == "TCR":
        t60_by_band, missing, invalid = GetT60PerOctaveBand_TCR(receiver_folder_id, gid)
        if t60_by_band is None:
            print(ui._("Could not read acoustic parameters for receiver '%s'.") % recv_name)
            return None
        if missing or invalid:
            if invalid:
                print(ui._(
                    "Receiver '%s': RT-30 is NaN/invalid in these IEC bands: %s. This "
                    "usually means too few particles reached the receiver in those "
                    "bands. Increase the particle count and re-run."
                ) % (recv_name, ", ".join("%d Hz" % hz for hz in invalid)))
            if missing:
                print(ui._(
                    "Receiver '%s': these IEC octave bands are not present in the "
                    "results: %s. Enable them in the calculation's frequency bands."
                ) % (recv_name, ", ".join("%d Hz" % hz for hz in missing)))
            print(ui._("STI needs all 7 bands (125 Hz - 8 kHz) and was NOT computed."))
            return None
        sti, mti_by_band = ComputeSTI_TCR(t60_by_band, reduction_factors)
    else:
        echogram_by_band, missing, decay_db = GetEchogramPerOctaveBand_SPPS(gid)
        if echogram_by_band is None:
            print(ui._("Could not read the echogram for receiver '%s'.") % recv_name)
            return None
        if missing:
            print(ui._(
                "Receiver '%s': these IEC octave bands are missing from the echogram: "
                "%s. STI needs all 7 bands (125 Hz - 8 kHz) and was NOT computed."
            ) % (recv_name, ", ".join("%d Hz" % hz for hz in missing)))
            return None
        sti, mti_by_band = ComputeSTI_SPPS(echogram_by_band, reduction_factors)
        worst_decay = min(decay_db.values()) if decay_db else None
        if worst_decay is not None and worst_decay < 20.0:
            truncation_warning = worst_decay

        # Second opinion from the Schroeder decay (see module notes).
        try:
            t60_sch, miss_sch, bad_sch = GetT60PerOctaveBand_TCR(
                receiver_folder_id, gid)
            if t60_sch and not miss_sch and not bad_sch:
                sti_schroeder, _mti_sch = ComputeSTI_TCR(
                    t60_sch, reduction_factors)
                t60_schroeder = t60_sch
        except Exception as exc:
            print(ui._("   (Schroeder/RT cross-check unavailable: %s)") % exc)

    return {
        "sti": sti,
        "mti_by_band": mti_by_band,
        "snr_db_by_band": snr_db_by_band,
        "speech_spl_by_band": speech_spl_by_band,
        "use_noise": use_noise,
        "truncation_warning": truncation_warning,
        "use_masking": use_masking,
        "reduction_factors": reduction_factors,
        "sti_schroeder": sti_schroeder,
        "t60_schroeder": t60_schroeder,
    }


def ComputeSTI(receiver_folder_id):
    receiver = ui.element(receiver_folder_id)
    recv_name = receiver.getinfos()["label"]

    mode = DetectCalculationMode(receiver_folder_id)
    if mode is None:
        print(ui._("Could not determine whether this belongs to SPPS or TCR results; aborting."))
        return

    gid = FindSoundLevelGabeId(receiver_folder_id)
    if gid == -1:
        print(ui._("Receiver '%s' has no 'Sound level' report; run a calculation first.") % recv_name)
        return

    # Background noise is specified the way IEC 60268-16 / ISO 9921 (and
    # do it: as background noise SPL per octave band. The SNR is then
    # DERIVED per band as speech level minus noise level, using the speech
    # level actually predicted at this receiver -- not guessed at.
    lbl_mode = ui._("Background noise")
    fields = {lbl_mode: list(NOISE_MODES),
              ui._(MASKING_LABEL): list(MASKING_MODES)}
    for hz in NOISE_ENTRY_BANDS_HZ:
        fields[_noise_label(hz)] = "%g" % DEFAULT_NOISE_SPL_DB

    res = ui.application.getuserinput(
        ui._("Compute STI for '%s' (%s branch)") % (recv_name, mode),
        ui._(
            "IEC 60268-16 full STI. Choose the background noise: none, an NC "
            "curve, or User Defined.\n\n"
            "IMPORTANT: the 'Noise SPL' fields below are used ONLY when "
            "'User Defined' is selected. Picking an NC curve does NOT rewrite "
            "them (this dialog cannot update its own fields) -- the NC levels "
            "are applied internally and the exact spectrum used is printed "
            "with the result and stored in the results file.\n\n"
            "The per-band SNR is computed from the speech level predicted at "
            "this receiver minus the background noise level."
        ),
        fields,
    )
    if not res[0]:
        return
    noise_mode = res[1].get(lbl_mode, NOISE_MODE_NONE)
    use_masking = res[1].get(ui._(MASKING_LABEL), MASKING_ON) == MASKING_ON

    noise_spl_by_band = {}
    noise_description = noise_mode
    if noise_mode == NOISE_MODE_USER:
        try:
            for hz in sm.OCTAVE_BANDS_HZ:
                noise_spl_by_band[hz] = float(res[1][_noise_label(hz)])
        except (ValueError, KeyError):
            print(ui._("Invalid background noise level."))
            return
    elif noise_mode != NOISE_MODE_NONE:
        nc_value = None
        for nc in sm.available_nc_values():
            if noise_mode == _nc_mode_label(nc):
                nc_value = nc
                break
        if nc_value is None:
            print(ui._("Unrecognised background noise selection: %s") % noise_mode)
            return
        noise_spl_by_band = sm.nc_curve_levels(nc_value)
        noise_description = u"NC %d curve" % nc_value

    _sti_result = sti_for_receiver(
        receiver_folder_id, gid, mode, recv_name, noise_spl_by_band,
        use_masking)
    if _sti_result is None:
        return
    sti = _sti_result["sti"]
    mti_by_band = _sti_result["mti_by_band"]
    snr_db_by_band = _sti_result["snr_db_by_band"]
    speech_spl_by_band = _sti_result["speech_spl_by_band"]
    use_noise = _sti_result["use_noise"]
    truncation_warning = _sti_result["truncation_warning"]
    use_masking = _sti_result["use_masking"]
    reduction_factors = _sti_result["reduction_factors"]
    sti_schroeder = _sti_result["sti_schroeder"]
    t60_schroeder = _sti_result["t60_schroeder"]

    rating_for_file, _pb = sm.iso9921_rating(sti)
    # The results spreadsheet only stores numbers in the value column, so
    # the verdict is carried in the ROW LABEL where it reads as words. The
    # numeric value stays meaningful (an ordinal) so it can still be
    # sorted or mapped.
    row_labels = ([
        "STI",
        "Rating (ISO 9921): %s" % rating_for_file.upper(),
        "Noise included: %s" % ("Yes" if use_noise else "No"),
        "Auditory masking applied: %s" % ("Yes" if use_masking else "No"),
        "Background noise: %s" % (noise_description if use_noise else "none"),
        "STI from Schroeder decay (RT-30)",
        "STI echogram minus Schroeder",
    ]
        + ["MTI %d Hz" % hz for hz in sm.OCTAVE_BANDS_HZ]
        + ["SNR %d Hz (dB)" % hz for hz in sm.OCTAVE_BANDS_HZ])
    # NB: pass plain str, NOT bytes. libsimpa's SWIG binding maps
    # std::string to Python 3 str; passing the .encode('cp1252') bytes that
    # the shipped Python 2-era sample scripts use raises
    # "TypeError: in method 'stringarray_append', argument 2 of type
    # 'std::vector< std::string >::value_type const &'" (confirmed live).
    gabewriter = Gabe_rw(2)
    labelcol = stringarray()
    for lbl in row_labels:
        labelcol.append(lbl)
    gabewriter.AppendStrCol(labelcol, "label")
    rating_order = ["bad", "poor", "fair", "good", "excellent"]
    datacol = floatarray()
    datacol.append(float(sti))
    datacol.append(float(rating_order.index(rating_for_file)))
    datacol.append(1.0 if use_noise else 0.0)
    datacol.append(1.0 if use_masking else 0.0)
    # NC number if an NC curve was used, else 0
    nc_used = 0.0
    for nc in sm.available_nc_values():
        if noise_mode == _nc_mode_label(nc):
            nc_used = float(nc)
            break
    datacol.append(nc_used)
    # order must match row_labels: the two Schroeder rows come next
    datacol.append(float("nan") if sti_schroeder is None else float(sti_schroeder))
    datacol.append(float("nan") if sti_schroeder is None
                   else float(sti - sti_schroeder))
    for hz in sm.OCTAVE_BANDS_HZ:
        datacol.append(float(mti_by_band[hz]))
    for hz in sm.OCTAVE_BANDS_HZ:
        datacol.append(float(snr_db_by_band[hz]))
    gabewriter.AppendFloatCol(datacol, recv_name)

    out_path = receiver_output_path(receiver_folder_id, "STI.gabe")
    gabewriter.Save(out_path)
    ui.application.sendevent(receiver, ui.idevent.IDEVENT_RELOAD_FOLDER)

    rating, pb_score = sm.iso9921_rating(sti)

    print(ui._("=== STI for '%s' (%s branch) ===") % (recv_name, mode))
    _n_per_source = CountPerSourceEchograms(receiver_folder_id)
    if _n_per_source:
        print(ui._(
            "  Using the COMBINED echogram (all sources summed). "
            "'Echogram per source' is enabled, so %d extra per-source "
            "echogram(s) also exist below this receiver; they are NOT used."
        ) % _n_per_source)
    print("STI = %.3f   ->  %s (ISO 9921)"
          % (sti, rating.upper()))
    print(ui._("   meaningful PB-word score: %s ; JND for STI is %.2f")
          % (pb_score, sm.STI_JND))
    if sti_schroeder is not None:
        delta = sti - sti_schroeder
        print(ui._(
            "   Cross-check from the SCHROEDER decay: STI = %.3f "
            "(%s, %+.3f vs the echogram result)")
            % (sti_schroeder, sm.iso9921_rating(sti_schroeder)[0].upper(), delta))
        print(ui._(
            "      That estimate uses I-Simpa's RT-30, fitted by linear "
            "regression to the Schroeder backward integration, with the "
            "analytical exponential MTF m(F)=1/sqrt(1+(2*pi*F*tau)^2), "
            "tau=T60/13.8. It therefore assumes a PURELY EXPONENTIAL decay."))
        if t60_schroeder:
            print("      " + "".join("%9d" % hz for hz in sm.OCTAVE_BANDS_HZ) + "   (Hz)")
            print("      " + "".join("%9.3f" % t60_schroeder[hz]
                                     for hz in sm.OCTAVE_BANDS_HZ) + "   (RT-30, s)")
        if abs(delta) <= sm.STI_JND:
            print(ui._(
                "      The two agree to within one JND, so this room's decay "
                "is close to a single exponential and both numbers stand."))
        elif delta > 0:
            print(ui._(
                "      The echogram result is HIGHER. Early/direct energy is "
                "helping intelligibility in a way a single decay slope cannot "
                "represent -- expected with several sources, or a receiver "
                "close to one. The echogram value is the better estimate; the "
                "Schroeder one is the conservative floor."))
        else:
            print(ui._(
                "      The echogram result is LOWER. Distinct echoes or a "
                "double-sloped decay are destroying modulation in a way RT "
                "alone cannot see. Trust the echogram value and inspect the "
                "echogram for late reflections."))
    if use_noise:
        print(ui._("   background noise: %s ; per-band speech / noise / SNR (dB):")
              % noise_description)
        print("      " + "".join("%9d" % hz for hz in sm.OCTAVE_BANDS_HZ) + "   (Hz)")
        print("      " + "".join("%9.1f" % speech_spl_by_band[hz] for hz in sm.OCTAVE_BANDS_HZ) + "   (speech)")
        print("      " + "".join("%9.1f" % noise_spl_by_band[hz] for hz in sm.OCTAVE_BANDS_HZ) + "   (noise)")
        print("      " + "".join("%9.1f" % snr_db_by_band[hz] for hz in sm.OCTAVE_BANDS_HZ) + "   (SNR)")
    else:
        print(ui._(
            "   NOISE-FREE result: reverberation only, no background noise. This is "
            "an upper bound -- any real background noise will lower it."))
    print("  " + "".join("%9d" % hz for hz in sm.OCTAVE_BANDS_HZ) + "   (Hz)")
    print("  " + "".join("%9.3f" % mti_by_band[hz] for hz in sm.OCTAVE_BANDS_HZ) + "   (MTI)")
    if truncation_warning is not None:
        print(ui._(
            "WARNING: the echogram only decays by %.1f dB before it ends, so the "
            "impulse response is truncated well above the noise floor. The MTF "
            "integral then misses late reverberant energy, which biases STI "
            "HIGH (too optimistic). Increase the calculation duration so the "
            "decay covers at least 20-30 dB before trusting this value."
        ) % truncation_warning)
    print(ui._(
        "   NOTE - amplified vs unamplified: STI itself does not know which "
        "this is; that is set by the SOURCE you defined. For an UNAMPLIFIED "
        "talker, set the source sound power to the ANSI S3.5 speech levels for "
        "the intended vocal effort (Normal / Raised / Loud / Shouted). For an "
        "AMPLIFIED system, use the loudspeaker sensitivity and directivity, and "
        "per IEC 60268-16 apply the male speech spectrum as the equaliser. With "
        "no background noise the two give IDENTICAL STI, because the modulation "
        "transfer function is level-independent."))
    print(ui._("Wrote: %s") % out_path)


class manager:
    def __init__(self):
        self.compute_sti_id = ui.application.register_event(self.OnComputeSTI)

    def getmenu(self, typeel, idel, menu):
        el = ui.element(idel)
        infos = el.getinfos()
        parentid = infos.get("parentid", -1)
        if parentid is None or parentid < 0:
            return False
        parent = ui.element(parentid)
        parent_infos = parent.getinfos()
        if parent_infos.get("name") == PUNCTUAL_RECEIVERS_FOLDER_NAME:
            menu.insert(0, ())
            menu.insert(0, (ui._("Compute STI"), self.compute_sti_id))
            return True
        return False

    def OnComputeSTI(self, idel):
        ComputeSTI(idel)


ui.application.register_menu_manager(ui.element_type.ELEMENT_TYPE_REPORT_FOLDER, manager())
