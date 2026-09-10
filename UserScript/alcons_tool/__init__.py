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
# %ALcons (Peutz articulation loss of consonants) post-processing tool.
#
# Adds a "Compute %ALcons" entry to the right-click menu of each INDIVIDUAL
# punctual receiver folder (e.g. "Receiver 1").
#
# ONE METHOD ONLY: the ENERGETIC form, which is what commercial
# tools use because it works for any number of sources with any
# directivities:
#
#       %ALcons = 8.9 * T60 * E_rev / (13.82 * E_dir)
#
# (Bies & Hansen, "Engineering Noise Control", eq 7.130, after Bistafa &
# Bradley 2000. The text recommends it precisely because the classical
# formula "is difficult to use when there are multiple sound sources such
# as in a typical sound reinforcement system".)
#
# WHY THERE IS NO Q, NO V, NO r AND NO SOURCE PICKER
# --------------------------------------------------
# The energetic form needs none of them. E_dir and E_rev are read from the
# SPPS echogram at this receiver, and SPPS has ALREADY applied each
# source's own directivity (balloon / XY / XZ / YZ plan / unidirectional /
# omnidirectional) when it launched that source's particles -- see
# sppsNantes.cpp, which multiplies each particle's energy by the balloon
# value for its launch direction. Whatever sources exist in the project,
# and whatever directivity each one is set to, is exactly what is used.
#
# Deliberately NOT offering a Q input also removes any way for a user to
# silently contradict the directivity actually configured on the sources.
#
# ON THE CONSTANT (please read before changing it)
# ------------------------------------------------
# A "%ALcons = 9 * (E_rev/E_dir) * T60" form circulates widely. It is the
# same equation but is missing the 13.82 divisor, and overestimates by
# ~13.8x. Three independent checks, all agreeing:
#
#  1. Substituting the steady-state relations E_rev/E_dir = 16*pi*r^2/(Q*R)
#     and Sabine A = 0.161V/T60 into eq 7.130 yields
#     201.1 * r^2 * T60^2 / (V*Q), i.e. it reproduces the classical
#     Peutz coefficient of 200 to within 0.5%. The 9x form yields 2810.
#  2. At the critical distance rH, E_dir = E_rev by definition, so both
#     forms must agree there. Classical eq 7.129 evaluated at rH gives
#     0.65*T60; eq 7.130 gives 0.644*T60 (-0.9%); the 9x form gives
#     9*T60 (+1285%).
#  3. 13.82 = ln(10^6), the T60 -> decay-constant conversion. The 9x form
#     is only correct if its "T60" is really tau = T60/13.82.
#
# If you nevertheless want the 9x scaling, set NUMERATOR = 9.0 and
# DENOMINATOR = 1.0 below -- it is a one-line change and the active
# formula is printed with every result so it can never be ambiguous.
#
# GRID LAYOUT (verified live via debug dump, not assumed):
#   'Acoustic parameters': row[0] = ['', ..., 'RT-30 (s)', 'EDT (s)', ...]
#                          row[i] = ['63 Hz', <values>]
#   'Sound level'        : row[0] = ['', '10.0 ms', ..., 'Total']
#                          row[i] = ['63 Hz', <one value per time step>]
#   i.e. ROWS are bands, COLUMNS are time steps / parameters.
from __future__ import print_function
import math
import re
import uictrl as ui

# The empirical (Becker) %ALcons cross-check needs an STI for this same
# receiver. It is computed by sti_tool, reusing the very same function
# sti_tool's own report uses, so the two can never disagree. If sti_tool is
# not installed the cross-check is simply skipped -- the primary,
# energy-based %ALcons does not depend on it.
try:
    import sti_tool as _sti_tool
    import sti_math as _sm
    import sti_noise as _sti_noise
except Exception:
    _sti_tool = None
    _sm = None
    _sti_noise = None
from libsimpa import *

TARGET_BAND_HZ = 2000.0  # Peutz's basis; text: "the 2000 Hz octave band"
PUNCTUAL_RECEIVERS_FOLDER_NAME = u"Punctual receivers"
BAND_MATCH_TOLERANCE_HZ = 1.0

# Speed of sound is derived from the project's own Environment temperature
# using I-Simpa's OWN formula, so arrival times agree exactly with how the
# engine propagated the particles:
#   CCalculsGenerauxThermodynamique::c_son(K) = 343.2 * sqrt(K / 293.15)
# with K = temperature(degC) + 273.15   (Celerite_du_son.cpp, and
# calculsPropagation.h: Kref = 293.15). At 20 degC this gives exactly
# 343.2 m/s.
CELERITY_REF_MS = 343.2
CELERITY_REF_K = 293.15
DEFAULT_TEMPERATURE_C = 20.0

# --- coefficients, and how the two published constants fit together ---
#
# Bistafa & Bradley (2000) ratio term, Bies & Hansen eq (7.130):
ALCONS_NUMERATOR = 8.9
ALCONS_DENOMINATOR = 13.82
# Peutz's ceiling: "without direct sound, ALcons is limited to 9 * T60".
ALCONS_CEILING_FACTOR = 9.0
#
# These are NOT competing conventions -- they are two parts of one formula:
#   E_rev/E_dir grows as (r/Dc)^2 (it is 1 at the critical distance Dc by
#   definition). The classical Peutz expression reaches the 9*T60 ceiling at
#   E_rev/E_dir = 13.85, i.e. r = 3.72*Dc, i.e. D/R = -11.4 dB -- which
#   matches the textbook statement that D/R = -11 dB at 3.5*Dc. A ratio-based
#   form C*(E_rev/E_dir)*T60 must therefore use C = 9/13.85 = 0.650 to reach
#   the ceiling in the right place, and Bistafa & Bradley's 8.9/13.82 = 0.644
#   IS that number.
#
# So the "9" seen in circulating versions of this equation is real, but it is
# the CEILING constant, not the multiplier on the ratio. Using it as the
# multiplier reaches the ceiling at Dc instead of 3.72*Dc and then exceeds it,
# which is why it returns impossible (>100%) values in reverberant rooms.
#
# Both are used here: the ratio term below, saturated at Peutz's ceiling.

# eq (7.131) direct-to-reverberant level rating, in 3 dB steps
# (> -3 dB excellent ... < -15 dB very poor)
DR_RATING_BANDS_HZ = [500.0, 1000.0, 2000.0]
DR_RATING_STEPS = [
    (-3.0, "excellent"),
    (-6.0, "very good"),
    (-9.0, "good"),
    (-12.0, "fair"),
    (-15.0, "poor"),
]
DR_RATING_WORST = "very poor"

METHOD_BISTAFA = u"Bistafa & Bradley (2000) - ratio, capped at Peutz 9*T60"
METHOD_PEUTZ_LIMIT = u"Peutz limit only - 9*T60 (no-direct-sound case)"

_BAND_RE = re.compile(r"(\d+(?:[.,]\d+)?)\s*(k?)Hz", re.IGNORECASE)


def formula_text():
    # NB: '%%ALcons' -- a literal percent sign must be escaped here or the
    # '%A' is parsed as a format specifier (ValueError).
    return "%%ALcons = min( %g * T60 * E_rev / (%g * E_dir) , %g * T60 )" % (
        ALCONS_NUMERATOR, ALCONS_DENOMINATOR, ALCONS_CEILING_FACTOR)


def rate_direct_to_reverberant(dr_db):
    """Bies & Hansen eq (7.131) rating scale, in 3 dB steps."""
    for threshold, label in DR_RATING_STEPS:
        if dr_db > threshold:
            return label
    return DR_RATING_WORST


def _distance(p1, p2):
    return math.sqrt(sum((p1[i] - p2[i]) ** 2 for i in range(3)))


def speed_of_sound_from_temperature(temperature_c):
    """I-Simpa's own celerity law (Celerite_du_son.cpp):
    c = 343.2 * sqrt((T + 273.15) / 293.15)"""
    return CELERITY_REF_MS * math.sqrt((temperature_c + 273.15) / CELERITY_REF_K)


def GetProjectTemperature():
    """Reads Temperature (degC) from the project's Environment element.
    Returns (temperature_c, was_found)."""
    try:
        rootscene = ui.element(ui.application.getrootscene())
        env_id = rootscene.getelementbytype(
            ui.element_type.ELEMENT_TYPE_SCENE_PROJET_ENVIRONNEMENTCONF)
        if env_id is not None and env_id != -1:
            env = ui.element(env_id)
            if env.hasproperty("temperature"):
                value = _as_finite_float(env.getdecimalconfig("temperature"))
                if value is not None:
                    return value, True
    except Exception:
        pass
    return DEFAULT_TEMPERATURE_C, False


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
    """'10.0 ms' -> 0.010. None for 'Total' / the empty corner cell."""
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


def _as_finite_float(value):
    """None for NaN / inf / non-numeric. Low particle counts legitimately
    produce NaN RT values in some bands and those must be caught."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f != f or f == float("inf") or f == float("-inf"):
        return None
    return f


def _db_to_linear_energy(value):
    """Converts an echogram cell (dB) to linear energy.

    -inf dB means literally zero energy in that time step, which I-Simpa
    emits for steps before the direct sound arrives (confirmed in a live
    dump: the first steps of every band read '-inf'). That is VALID data
    meaning zero, NOT a failure, so it maps to 0.0 and the time step is
    kept -- dropping those rows would leave a ragged time axis. NaN or a
    non-numeric cell is a genuine failure and returns None."""
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if f != f:            # NaN -> invalid
        return None
    if f == float("-inf"):
        return 0.0
    if f == float("inf"):
        return None
    return 10 ** (f / 10.0)


def receiver_output_path(receiver_folder_id, filename):
    """Folder path helper. buildfullpath() lives on e_file, NOT on the
    plain element wrapper -- calling it on element raises AttributeError
    (confirmed live)."""
    return ui.e_file(receiver_folder_id).buildfullpath() + filename


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


SOURCE_ENABLE_PROPERTY = "enable"   # labelled "Active source" in the GUI


def SourceIsActive(source_id):
    """Whether this source is ticked "Active source".

    The flag is a BOOL on the source's PROPERTIES child, not on the source
    element itself, and Element::GetBoolConfig only looks at DIRECT
    children (element.cpp:1300-1314) -- so the children are walked.
    hasproperty() is used to tell "flag absent" from "flag is false",
    because GetBoolConfig returns false for BOTH (element.cpp:1312).

    Unknown means active, which preserves the old behaviour for any
    project whose sources carry no flag."""
    try:
        src = ui.element(source_id)
    except Exception:
        return True
    try:
        if src.hasproperty(SOURCE_ENABLE_PROPERTY):
            return bool(src.getboolconfig(SOURCE_ENABLE_PROPERTY))
    except Exception:
        pass
    try:
        children = src.childs()
    except Exception:
        return True
    for child in children:
        try:
            el = ui.element(child[0])
            if el.hasproperty(SOURCE_ENABLE_PROPERTY):
                return bool(el.getboolconfig(SOURCE_ENABLE_PROPERTY))
        except Exception:
            continue
    return True


def GetSceneSourcePositions(include_inactive=False):
    """Positions of the ACTIVE sources. Directivity does not need to be
    read here -- SPPS already baked it into the echogram.

    Inactive sources are excluded because they are not in the calculation
    at all (see this module's notes on SourceIsActive); giving one a direct
    arrival steals energy from E_rev into E_dir."""
    rootscene = ui.element(ui.application.getrootscene())
    src_ids = rootscene.getallelementbytype(ui.element_type.ELEMENT_TYPE_SCENE_SOURCES_SOURCE)
    positions = {}
    for sid in src_ids:
        src = ui.element(sid)
        name = src.getinfos()["name"]
        if not include_inactive and not SourceIsActive(sid):
            continue
        pos = src.getpositionconfig("pos_source")
        positions[name] = (pos[0], pos[1], pos[2])
    return positions


def GetInactiveSceneSourceNames():
    """Names of sources present in the scene but NOT in the calculation."""
    rootscene = ui.element(ui.application.getrootscene())
    try:
        src_ids = rootscene.getallelementbytype(
            ui.element_type.ELEMENT_TYPE_SCENE_SOURCES_SOURCE)
    except Exception:
        return []
    names = []
    for sid in src_ids:
        if SourceIsActive(sid):
            continue
        try:
            names.append(ui.element(sid).getinfos()["name"])
        except Exception:
            continue
    return names


def GetSceneReceiverPosition(receiver_name):
    rootscene = ui.element(ui.application.getrootscene())
    recp_group_id = rootscene.getelementbytype(ui.element_type.ELEMENT_TYPE_SCENE_RECEPTEURSP)
    if recp_group_id == -1:
        return None
    recp_ids = ui.element(recp_group_id).getallelementbytype(
        ui.element_type.ELEMENT_TYPE_SCENE_RECEPTEURSP_RECEPTEUR
    )
    for rid in recp_ids:
        recp = ui.element(rid)
        if recp.getinfos()["name"] == receiver_name:
            pos = recp.getpositionconfig("pos_recepteur")
            return (pos[0], pos[1], pos[2])
    return None


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


def GetT60Table(receiver_folder_id, sound_level_gabe_id):
    """Triggers I-Simpa's own acoustic-parameter calculation (so RT does
    not have to be computed by hand first), then returns
    (valid_t60_by_band_hz, invalid_band_list)."""
    ui.application.sendevent(
        ui.element(sound_level_gabe_id),
        ui.idevent.IDEVENT_RECP_COMPUTE_ACOUSTIC_PARAMETERS,
        {"TR": "30"},
    )
    receiver = ui.element(receiver_folder_id)
    params_id = receiver.getelementbylibelle("Acoustic parameters")
    if params_id == -1:
        return None, None
    grid = ui.application.getdataarray(ui.element(params_id))
    if not grid or len(grid) < 2:
        return None, None
    header = grid[0]
    idcol = _find_rt_column(header, 30)
    if idcol == -1:
        print(ui._("Could not find an RT-30 column. Header was: %r") % (header,))
        return None, None
    valid = {}
    invalid = []
    for row in grid[1:]:
        hz = _band_label_to_hz(row[0])
        if hz is None or idcol >= len(row):
            continue
        t60 = _as_finite_float(row[idcol])
        if t60 is None or t60 <= 0:
            invalid.append(hz)
        else:
            valid[hz] = t60
    return valid, invalid


def GetEchogramForBand(sound_level_gabe_id, target_hz):
    """Returns (times_s, energies_linear, matched_band_hz).

    Grid values are dB; converted to linear via 10**(dB/10) because the
    direct/reverberant split is an ENERGY ratio."""
    grid = ui.application.getdataarray(ui.element(sound_level_gabe_id))
    if not grid or len(grid) < 2:
        return None, None, None
    header = grid[0]

    time_cols = []
    for idx in range(1, len(header)):
        t = _time_label_to_seconds(header[idx])
        if t is not None:
            time_cols.append((t, idx))
    if len(time_cols) < 2:
        return None, None, None
    time_cols.sort(key=lambda pair: pair[0])

    best_row = best_hz = best_diff = None
    for row in grid[1:]:
        hz = _band_label_to_hz(row[0])
        if hz is None:
            continue
        diff = abs(hz - target_hz)
        if best_diff is None or diff < best_diff:
            best_diff, best_row, best_hz = diff, row, hz
    if best_row is None or best_diff > BAND_MATCH_TOLERANCE_HZ:
        return None, None, None

    times, energies = [], []
    for t, idx in time_cols:
        if idx >= len(best_row):
            continue
        e = _db_to_linear_energy(best_row[idx])
        if e is None:
            continue
        times.append(t)
        energies.append(e)
    if len(times) < 2 or sum(energies) <= 0:
        return None, None, None
    return times, energies, best_hz


def SplitDirectReverberant(times, energies, arrival_times_s):
    """Splits echogram energy into direct and reverberant parts.

    Time-step labels are the END of each step (see
    BaseReportManager::InitHeaderArrays), so step k covers
    (times[k]-dt, times[k]]. Each source's direct sound lands in the step
    containing its geometric arrival time r/c. Every step containing any
    source's direct arrival counts toward E_dir (a step shared by two
    sources is counted once); everything else is E_rev.

    Returns (E_dir, E_rev, direct_step_times, dt).
    """
    dt = times[0] if len(times) == 1 else min(
        times[i + 1] - times[i] for i in range(len(times) - 1)
    )
    direct_indices = set()
    for t_arr in arrival_times_s:
        for k, t_end in enumerate(times):
            if t_arr <= t_end + 1e-12:
                # The step containing the geometric arrival can still hold
                # zero energy (no particle has landed yet -- I-Simpa writes
                # -inf dB there). The direct sound is then the first step
                # that actually carries energy, so walk forward to it
                # rather than reporting E_dir = 0.
                while k < len(energies) and energies[k] <= 0.0:
                    k += 1
                if k < len(energies):
                    direct_indices.add(k)
                break
    e_dir = sum(energies[k] for k in direct_indices)
    e_rev = sum(e for k, e in enumerate(energies) if k not in direct_indices)
    return e_dir, e_rev, sorted(times[k] for k in direct_indices), dt


def AlconsRatioTerm(t60, e_dir, e_rev):
    """Bies & Hansen eq (7.130) / Bistafa & Bradley (2000), UNCAPPED."""
    if e_dir <= 0:
        return None
    return (ALCONS_NUMERATOR * t60 * e_rev) / (ALCONS_DENOMINATOR * e_dir)


def AlconsPeutzCeiling(t60):
    """Peutz: without direct sound, %ALcons is limited to 9 * T60."""
    return ALCONS_CEILING_FACTOR * t60


def AlconsEnergetic(t60, e_dir, e_rev):
    """The full Bistafa & Bradley result: the ratio term, saturated at
    Peutz's no-direct-sound ceiling. Returns (alcons, is_saturated)."""
    ratio_term = AlconsRatioTerm(t60, e_dir, e_rev)
    if ratio_term is None:
        return None, False
    ceiling = AlconsPeutzCeiling(t60)
    if ratio_term >= ceiling:
        return ceiling, True
    return ratio_term, False


def _alcons_dialog_fields(lbl_method):
    """Algorithm choice, plus the background-noise spectrum when the STI
    cross-check is available to use it."""
    fields = {lbl_method: [METHOD_BISTAFA, METHOD_PEUTZ_LIMIT]}
    if _sti_noise is not None:
        fields.update(_sti_noise.dialog_fields())
    return fields


def _becker_cross_check(receiver_folder_id, gid, recv_name,
                        noise_spl_by_band, alcons_energy):
    """Secondary %ALcons, derived from this receiver's STI.

    Returns (sti, alcons_percent, trusted, ratio_to_energy_result), with
    None entries when it could not be computed."""
    if _sti_tool is None or _sm is None:
        return None, None, None, None
    try:
        mode = _sti_tool.DetectCalculationMode(receiver_folder_id)
        if mode is None:
            return None, None, None, None
        result = _sti_tool.sti_for_receiver(
            receiver_folder_id, gid, mode, recv_name, noise_spl_by_band)
    except Exception as exc:
        print(ui._("  (STI cross-check unavailable: %s)") % exc)
        return None, None, None, None
    if not result:
        return None, None, None, None
    sti = result["sti"]
    alcons_sti, trusted = _sm.alcons_from_sti(sti)
    divergence = None
    if alcons_energy and alcons_energy > 0:
        divergence = alcons_sti / alcons_energy
    return sti, alcons_sti, trusted, divergence


# A factor of two between two estimates of the same quantity is well beyond
# the spread of the published fits, so it points at a real modelling issue
# (non-diffuse field, or noise dominating) rather than formula scatter.
DIVERGENCE_FACTOR = 2.0


def _report_cross_check(sti, alcons_sti, trusted, divergence,
                        alcons_energy, noise_description):
    if sti is None:
        print(ui._(
            "  Cross-check: not available (sti_tool not installed, or STI "
            "could not be computed for this receiver)."))
        return
    print(ui._("  SECONDARY cross-check, %ALcons derived from STI:"))
    print("      STI = %.3f  ->  %%ALcons = %.2f %%   (background noise: %s)"
          % (sti, alcons_sti, noise_description))
    print(ui._(
        "      Empirical conversion 170.54*exp(-5.419*STI), attributed to "
        "Farrel Becker. NOT part of IEC 60268-16, and its own publisher "
        "notes it should be treated with some doubt. Reported only as a "
        "cross-check; the primary result above is the energy-based one."))
    if not trusted:
        print(ui._(
            "      NOTE: below STI 0.20 this fit is extrapolation and the "
            "value is capped at 100 %; treat it as indicative only."))
    if divergence is not None:
        print("      Becker / energy = %.2f x" % divergence)
        if divergence >= DIVERGENCE_FACTOR or divergence <= 1.0 / DIVERGENCE_FACTOR:
            print(ui._(
                "      *** The two estimates differ by more than a factor of "
                "%g. They measure the same thing by different routes, so a gap "
                "this large usually means the field is not diffuse (strong "
                "direct sound or distinct echoes), or that background noise "
                "rather than reverberation is limiting intelligibility -- the "
                "energy method cannot see noise at all. Trust the STI-derived "
                "value when noise dominates."
            ) % DIVERGENCE_FACTOR)


def ComputeALcons(receiver_folder_id):
    receiver = ui.element(receiver_folder_id)
    recv_name = receiver.getinfos()["label"]

    positions_src = GetSceneSourcePositions()
    if not positions_src:
        if GetInactiveSceneSourceNames():
            print(ui._(
                "Every sound source in this scene is unticked ('Active "
                "source' is off), so nothing was radiated and %ALcons cannot "
                "be computed. Activate at least one source and re-run."))
        else:
            print(ui._("No sound sources found in the scene."))
        return
    recv_pos = GetSceneReceiverPosition(recv_name)
    if recv_pos is None:
        print(ui._("No matching scene position found for receiver '%s'.") % recv_name)
        return

    gid = FindSoundLevelGabeId(receiver_folder_id)
    if gid == -1:
        print(ui._("Receiver '%s' has no 'Sound level' report; run a calculation first.") % recv_name)
        return

    # The only choice offered is WHICH PUBLISHED ALGORITHM to report. No
    # physical quantity is asked for: temperature comes from the project's
    # Environment element and everything else from the results, so nothing
    # a user sets here can contradict the model.
    lbl_method = ui._("Algorithm")
    res = ui.application.getuserinput(
        ui._("Compute %%ALcons for '%s'") % recv_name,
        ui._(
            "Direct and reverberant energy are taken from this receiver's "
            "echogram, so every source and its own directivity is already "
            "accounted for. Both algorithms below use Peutz's 9*T60 "
            "no-direct-sound limit; they differ in whether the "
            "direct/reverberant ratio is used as well.\n\n"
            "The energy method is BLIND TO BACKGROUND NOISE -- it only knows "
            "the direct/reverberant split. The background noise chosen below "
            "is therefore NOT used for the primary result; it feeds the "
            "secondary, STI-derived cross-check, which does account for it."
        ),
        _alcons_dialog_fields(lbl_method),
    )
    if not res[0]:
        return
    method = res[1].get(lbl_method, METHOD_BISTAFA)

    noise_spl_by_band, noise_description = {}, u"none"
    if _sti_noise is not None:
        ok, noise_spl_by_band, noise_description = _sti_noise.resolve(res[1])
        if not ok:
            return

    temperature_c, temp_found = GetProjectTemperature()
    speed_of_sound = speed_of_sound_from_temperature(temperature_c)

    t60_by_band, invalid_bands = GetT60Table(receiver_folder_id, gid)
    if t60_by_band is None:
        print(ui._("Could not read acoustic parameters for receiver '%s'.") % recv_name)
        return

    t60 = None
    for hz, value in t60_by_band.items():
        if abs(hz - TARGET_BAND_HZ) <= BAND_MATCH_TOLERANCE_HZ:
            t60 = value
            break
    if t60 is None:
        if any(abs(hz - TARGET_BAND_HZ) <= BAND_MATCH_TOLERANCE_HZ for hz in (invalid_bands or [])):
            print(ui._(
                "Receiver '%s': the 2000 Hz RT-30 is NaN/invalid, which normally means "
                "too few particles reached this receiver in that band. Increase the "
                "particle count and re-run. %%ALcons was NOT computed."
            ) % recv_name)
        else:
            print(ui._(
                "Receiver '%s': no 2000 Hz band in the results. %%ALcons uses the 2 kHz "
                "octave band and will not substitute another band."
            ) % recv_name)
        if t60_by_band:
            print(ui._("Bands with a valid RT-30: %s")
                  % ", ".join("%g Hz" % hz for hz in sorted(t60_by_band.keys())))
        return

    times, energies, band_hz = GetEchogramForBand(gid, TARGET_BAND_HZ)
    if times is None:
        print(ui._("Could not read a 2000 Hz echogram for receiver '%s'.") % recv_name)
        return

    arrivals = {}
    for name, pos in positions_src.items():
        r_i = _distance(pos, recv_pos)
        arrivals[name] = (r_i, r_i / speed_of_sound)

    e_dir, e_rev, direct_times, dt = SplitDirectReverberant(
        times, energies, [t for (_, t) in arrivals.values()])
    if e_dir <= 0:
        print(ui._("Could not identify any direct-sound arrival inside the echogram."))
        return

    alcons_bistafa, saturated = AlconsEnergetic(t60, e_dir, e_rev)
    alcons_ratio_only = AlconsRatioTerm(t60, e_dir, e_rev)
    alcons_ceiling = AlconsPeutzCeiling(t60)
    alcons = alcons_ceiling if method == METHOD_PEUTZ_LIMIT else alcons_bistafa

    ratio = e_rev / e_dir
    d_over_r_db = 10 * math.log10(e_dir / e_rev) if e_rev > 0 else float("inf")

    # eq (7.131): direct-to-reverberant level for 500 Hz, 1 kHz and 2 kHz
    dr_by_band = {}
    for band_target in DR_RATING_BANDS_HZ:
        b_times, b_energies, _b_hz = GetEchogramForBand(gid, band_target)
        if b_times is None:
            continue
        b_dir, b_rev, _bt, _bdt = SplitDirectReverberant(
            b_times, b_energies, [t for (_, t) in arrivals.values()])
        if b_dir > 0 and b_rev > 0:
            dr_by_band[band_target] = 10 * math.log10(b_dir / b_rev)

    sti_x, alcons_x, alcons_x_trusted, divergence = _becker_cross_check(
        receiver_folder_id, gid, recv_name, noise_spl_by_band, alcons)

    # Trimmed table. Three rows were removed as redundant rather than
    # wrong:
    #   "Direct-to-reverberant 2kHz (dB)" -- an exact duplicate of the
    #       "D/R 2000 Hz (dB)" row appended below;
    #   "%ALcons Peutz limit 9*T60 (%)"   -- just 9 * the T60 row;
    #   "%ALcons ratio term, uncapped (%)"-- equal to %ALcons unless the
    #       result saturated, which the saturation flag already states.
    #
    # ORDER: the STI-derived value is listed FIRST because for a
    # multi-source system it is the one to quote. The energy-based number
    # keeps its full name but is explicitly labelled a diagnostic: its
    # direct/reverberant split degrades as (number of sources x time step)
    # grows, and with many sources it reads optimistically however fine
    # the time step.
    row_labels = [
        "STI (same receiver)",
        "%ALcons from STI, Becker (%) - REPORT THIS",
        "%ALcons energy-based (%) - diagnostic only",
        "Agreement ratio (Becker / energy)",
        "Saturated at Peutz limit (1=yes)",
        "T60 2kHz (s)",
        "E_rev / E_dir",
        "Number of ACTIVE sources",
        "Echogram time step (ms)",
    ]
    values = [float("nan") if sti_x is None else sti_x,
              float("nan") if alcons_x is None else alcons_x,
              alcons,
              float("nan") if divergence is None else divergence,
              1.0 if saturated else 0.0, t60, ratio,
              float(len(positions_src)), dt * 1000.0]
    for band_target in DR_RATING_BANDS_HZ:
        row_labels.append("D/R %d Hz (dB)" % int(band_target))
        values.append(dr_by_band.get(band_target, float("nan")))

    gabewriter = Gabe_rw(2)
    labelcol = stringarray()
    for lbl in row_labels:
        labelcol.append(lbl)
    gabewriter.AppendStrCol(labelcol, "label")
    datacol = floatarray()
    for v in values:
        datacol.append(float(v))
    gabewriter.AppendFloatCol(datacol, recv_name)

    out_path = receiver_output_path(receiver_folder_id, "ALcons.gabe")
    gabewriter.Save(out_path)
    ui.application.sendevent(receiver, ui.idevent.IDEVENT_RELOAD_FOLDER)

    print(ui._("=== %%ALcons for '%s' ===") % recv_name)
    _n_per_source = CountPerSourceEchograms(receiver_folder_id)
    if _n_per_source:
        print(ui._(
            "  Using the COMBINED echogram (all sources summed). "
            "'Echogram per source' is enabled, so %d extra per-source "
            "echogram(s) also exist below this receiver; they are NOT used."
        ) % _n_per_source)
    print("  " + ui._("Algorithm: %s") % method)
    print("  " + formula_text())
    print("  %%ALcons = %.2f %%   (T60 2kHz = %.3f s, E_rev/E_dir = %.2f)"
          % (alcons, t60, ratio))
    print("      ratio term (uncapped) = %.2f %% ; Peutz limit 9*T60 = %.2f %%"
          % (alcons_ratio_only, alcons_ceiling))
    if saturated and method == METHOD_BISTAFA:
        print(ui._(
            "      -> saturated: the ratio term exceeds Peutz's no-direct-sound "
            "limit, so the limit is reported. Expected in a very reverberant room."))
    if dr_by_band:
        print(ui._("  Direct-to-reverberant level, eq (7.131):"))
        for band_target in DR_RATING_BANDS_HZ:
            if band_target in dr_by_band:
                value = dr_by_band[band_target]
                print("      %5d Hz : %+6.1f dB   %s"
                      % (int(band_target), value, rate_direct_to_reverberant(value)))
    if temp_found:
        print(ui._("  c = %.1f m/s from the project Environment temperature of %.1f degC")
              % (speed_of_sound, temperature_c))
    else:
        print(ui._(
            "  WARNING: could not read the project Environment temperature; assumed "
            "%.1f degC -> c = %.1f m/s. Check the Environment element."
        ) % (temperature_c, speed_of_sound))
    _inactive = GetInactiveSceneSourceNames()
    print(ui._(
        "  Sources: %d ACTIVE of %d in the scene (only active ones get a "
        "direct-sound step; inactive ones are not in the calculation at all)."
    ) % (len(positions_src), len(positions_src) + len(_inactive)))
    if _inactive:
        print(ui._("      excluded as inactive: %s")
              % ", ".join(sorted(_inactive)))
    else:
        print(ui._(
            "      none were detected as inactive -- if you DID untick some, "
            "run 'DEBUG: Dump scene sources' and send the output, because "
            "the 'Active source' flag is then not where this tool looks."))
    print(ui._("  %d source(s), all directivities as configured in the project:")
          % len(positions_src))
    for name in sorted(arrivals.keys()):
        r_i, t_i = arrivals[name]
        print("      %-24s r = %6.2f m   direct arrival = %6.1f ms" % (name, r_i, t_i * 1000.0))
    print(ui._("  Echogram steps counted as direct: %s")
          % ", ".join("%.1f ms" % (t * 1000.0) for t in direct_times))

    shortest_path = min(r for (r, _) in arrivals.values())
    path_window = speed_of_sound * dt
    if path_window > 0.5 * shortest_path:
        print(ui._(
            "  WARNING: the echogram time step is %.1f ms, i.e. a %.1f m path window, "
            "which is large compared with the %.1f m direct path. Early reflections "
            "fall inside the direct step, inflating E_dir and making %%ALcons look "
            "BETTER than it is. Use a finer time step for a trustworthy split."
        ) % (dt * 1000.0, path_window, shortest_path))
    if invalid_bands:
        print(ui._("  Note: RT-30 is NaN/invalid in these bands (unused here): %s")
              % ", ".join("%g Hz" % hz for hz in sorted(invalid_bands)))
    _report_cross_check(sti_x, alcons_x, alcons_x_trusted, divergence,
                        alcons, noise_description)
    print(ui._("  Wrote: %s") % out_path)


class manager:
    def __init__(self):
        self.compute_alcons_id = ui.application.register_event(self.OnComputeALcons)

    def getmenu(self, typeel, idel, menu):
        el = ui.element(idel)
        infos = el.getinfos()
        parentid = infos.get("parentid", -1)
        if parentid is None or parentid < 0:
            return False
        parent = ui.element(parentid)
        if parent.getinfos().get("name") == PUNCTUAL_RECEIVERS_FOLDER_NAME:
            menu.insert(0, ())
            menu.insert(0, (ui._("Compute %ALcons"), self.compute_alcons_id))
            return True
        return False

    def OnComputeALcons(self, idel):
        ComputeALcons(idel)


ui.application.register_menu_manager(ui.element_type.ELEMENT_TYPE_REPORT_FOLDER, manager())
