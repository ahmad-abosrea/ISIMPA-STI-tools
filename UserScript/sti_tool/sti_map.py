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
# STI MAP on a surface receiver.
#
# Adds "Compute STI Map" to the right-click menu of a surface-receiver
# results folder. Computes STI per surface cell and writes a new
# surface-receiver result file that I-Simpa renders as a colour map.
#
# STRUCTURE (taken from a live dump, not assumed):
#
#   Surface receiver                     ELEMENT_TYPE_REPORT_FOLDER
#     +-- 63 Hz    -> rs_cut.csbin       one file per octave band
#     +-- 125 Hz   -> rs_cut.csbin
#     ...
#     +-- 16000 Hz -> rs_cut.csbin
#     +-- Global   -> rs_cut.csbin
#
#   Each file: record_type 'SPL_STANDART', 200 time steps of 0.01 s,
#   810 nodes, 1 surface receiver ('Plane Receiver', xmlid 666), 1508
#   faces. Per face: ~198 records, each a (time step index, energy) pair.
#
# TWO THINGS THAT MATTER, both confirmed from the dump:
#
#   1. The records are SPARSE and do not start at step 0 -- the first is
#      typically step=2. So each record's time MUST come from
#      GetFaceTimeStep(), never from the record index.
#
#   2. The energies are LINEAR (values like 3.99e-07), NOT dB. This is the
#      opposite of the punctual-receiver 'Sound level' grid, which is in
#      dB and needs 10**(dB/10). Here the values feed the MTF integral
#      directly -- applying a dB conversion would be badly wrong.
#
# OUTPUT RECORD TYPE: rsurf supports only "SPL_STANDART", "SPL_GAIN", "TR"
# and "EDT". There is no STI type. SPL types are converted to dB for
# display, which would destroy a 0..1 scale, so the map is written as "TR"
# (a time-like quantity, displayed as the raw number). The colours and
# values are then correct but the legend will read seconds rather than
# STI. Changing OUTPUT_RECORD_TYPE below is a one-line change if a
# different type renders better.
from __future__ import print_function
import math
import os
import sys

import uictrl as ui
import libsimpa as ls

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)
import sti_math as sm
import sti_noise

SURFACE_FOLDER_NAMES = (u"Surface receiver", u"Surface receivers")
OUTPUT_RECORD_TYPE = "TR"
OUTPUT_FILENAME = "STI_map.csbin"
PROGRESS_EVERY = 250

# rsurf_data::Make() only accepts "SPL_STANDART" / "SPL_GAIN" / "TR" / "EDT",
# so the FILE has to be written as one of those (TR is used: it is a
# time-like type, displayed raw, with no dB conversion applied).
#
# But the GUI does not read the record type from the file when it draws --
# it reads the 'rstype' INTEGER property on the tree element:
#     RecordRecepteurSurfType = elConf->GetIntegerConfig("rstype");
# and RECEPTEURS_RECORD_TYPE_STI is enum ordinal 9 (SPL_STANDART=0,
# SPL_GAIN=1, TR=2, EDT=3, PRESSURE=4, CLARITY=5, DEFINITION=6, TS=7,
# ST=8, STI=9).
#
# Overwriting that property to 9 after the element appears makes the GUI
# treat the map as a true STI map:
#   * legend unit becomes blank instead of "s"   (Recepteurs_surfacique.cpp)
#   * the load menu entry reads "Speech Transmision Index (STI)"
#   * auto-scaling stays on the RAW values, because the dB conversion is
#     applied only for SPL_STANDART / SPL_GAIN
RSTYPE_PROPERTY = "rstype"
RSTYPE_STI = 9


def _band_folder_map(folder_id):
    """{octave_hz: child folder element id} for the band subfolders."""
    folder = ui.element(folder_id)
    bands = {}
    try:
        children = folder.childs()
    except Exception:
        return bands
    for child in children:
        try:
            infos = ui.element(child[0]).getinfos()
        except Exception:
            continue
        name = infos.get("name") or u""
        hz = _label_to_hz(name)
        if hz is not None:
            bands[hz] = child[0]
    return bands


def _label_to_hz(label):
    if not isinstance(label, (str, bytes)):
        return None
    if isinstance(label, bytes):
        try:
            label = label.decode("cp1252")
        except Exception:
            return None
    text = label.strip()
    if not text.lower().endswith("hz"):
        return None
    try:
        return int(round(float(text[:-2].strip().replace(",", "."))))
    except ValueError:
        return None


def _find_csbin(folder_element_id):
    """The surface-receiver result file inside a band folder."""
    try:
        base = ui.e_file(folder_element_id).buildfullpath()
    except Exception:
        return None
    folder = base if os.path.isdir(base) else os.path.dirname(base)
    if not os.path.isdir(folder):
        return None
    for entry in sorted(os.listdir(folder)):
        if entry.lower().endswith(".csbin"):
            return os.path.join(folder, entry)
    return None


def _load_band(path):
    data = ls.rsurf_data()
    if not ls.rsurf_io.Load(path, data):
        return None
    return data


def _modulation_tables(n_steps, timestep):
    """Precomputed exp(-j2*pi*F*t) split into cos/sin, indexed
    [modulation frequency][time step].

    Built once instead of per face: with ~1500 faces x 7 bands x 14
    modulation frequencies x ~198 records this is the difference between
    seconds and many minutes."""
    cos_t, sin_t = [], []
    for f_mod in sm.MODULATION_FREQS_HZ:
        c_row, s_row = [], []
        for step in range(n_steps + 1):
            angle = 2.0 * math.pi * f_mod * step * timestep
            c_row.append(math.cos(angle))
            s_row.append(math.sin(angle))
        cos_t.append(c_row)
        sin_t.append(s_row)
    return cos_t, sin_t


def _face_total_energy(data, rs_index, face_index):
    """Summed energy for one face in one band, or None if it is empty."""
    n_rec = data.GetFaceRecordCount(rs_index, face_index)
    if n_rec <= 0:
        return None
    total = 0.0
    for rec in range(n_rec):
        energy = data.GetFaceEnergy(rs_index, face_index, rec)
        if energy > 0:
            total += energy
    return total if total > 0 else None


def _face_mti(data, rs_index, face_index, cos_t, sin_t, snr_db,
              noise_spl_db=None, reduction_factor=None):
    """MTI for one face in one band, straight from its sparse records.

    If noise_spl_db is given, the SNR is derived PER CELL from that cell's
    own speech level (its summed energy) rather than using a single value
    for the whole map -- the speech level varies across the surface, so a
    shared SNR would be wrong everywhere except one spot."""
    n_rec = data.GetFaceRecordCount(rs_index, face_index)
    if n_rec <= 0:
        return None

    steps, energies = [], []
    total = 0.0
    for rec in range(n_rec):
        energy = data.GetFaceEnergy(rs_index, face_index, rec)
        if energy <= 0:
            continue
        # time comes from the record's OWN step index; records are sparse
        steps.append(data.GetFaceTimeStep(rs_index, face_index, rec))
        energies.append(energy)
        total += energy
    if total <= 0 or len(steps) < 2:
        return None

    if reduction_factor is None:
        if noise_spl_db is not None:
            speech_spl = sti_noise.spl_from_energy(total)
            if speech_spl is None:
                return None
            snr_db = sm.snr_from_levels(speech_spl, noise_spl_db)
        reduction_factor = 1.0 / (1.0 + 10 ** (-snr_db / 10.0))

    n_table = len(cos_t[0]) - 1
    ti_values = []
    for fi in range(len(sm.MODULATION_FREQS_HZ)):
        c_row = cos_t[fi]
        s_row = sin_t[fi]
        re_sum = 0.0
        im_sum = 0.0
        for step, energy in zip(steps, energies):
            if step > n_table:
                continue
            re_sum += energy * c_row[step]
            im_sum -= energy * s_row[step]
        m = math.sqrt(re_sum * re_sum + im_sum * im_sum) / total
        if m > 1.0:
            m = 1.0
        ti_values.append(sm.transmission_index(m * reduction_factor))
    return sm.mti_for_band(ti_values)


def ComputeSTIMap(folder_id, noise_spl_by_band=None, noise_description=None,
                  use_masking=True):
    folder_infos = ui.element(folder_id).getinfos()
    print(ui._("=== STI map for %r ===") % folder_infos.get("name"))

    if noise_spl_by_band is None:
        noise_spl_by_band = {}
    if noise_description is None:
        noise_description = sti_noise.MODE_NONE
    use_noise = bool(noise_spl_by_band)
    snr_db_by_band = {hz: sti_noise.NOISE_FREE_SNR_DB for hz in sm.OCTAVE_BANDS_HZ}

    bands = _band_folder_map(folder_id)
    missing = [hz for hz in sm.OCTAVE_BANDS_HZ if hz not in bands]
    if missing:
        print(ui._("Missing octave band folder(s): %s. STI needs 125 Hz - 8 kHz.")
              % ", ".join("%d Hz" % hz for hz in missing))
        return

    loaded = {}
    for hz in sm.OCTAVE_BANDS_HZ:
        path = _find_csbin(bands[hz])
        if path is None:
            print(ui._("No .csbin result file found for the %d Hz band.") % hz)
            return
        data = _load_band(path)
        if data is None:
            print(ui._("Could not read %s") % path)
            return
        loaded[hz] = data
        print("   %5d Hz : %s" % (hz, os.path.basename(path)))

    reference = loaded[sm.OCTAVE_BANDS_HZ[0]]
    rs_count = reference.GetRsCount()
    n_steps = reference.GetTimeStepCount()
    timestep = reference.GetTimeStep()
    node_count = reference.GetNodesCount()

    # every band must describe the same geometry, or faces would not line up
    for hz, data in loaded.items():
        if (data.GetRsCount() != rs_count
                or data.GetNodesCount() != node_count):
            print(ui._("Band %d Hz has a different geometry; cannot combine bands.") % hz)
            return

    print(ui._("   %d surface receiver(s), %d nodes, %d time steps of %g s")
          % (rs_count, node_count, n_steps, timestep))

    cos_t, sin_t = _modulation_tables(n_steps, timestep)

    out = ls.rsurf_data()
    out.Make(node_count, rs_count, 1, timestep, OUTPUT_RECORD_TYPE)
    for node in range(node_count):
        pos = reference.GetNodePositionValue(node)
        out.SetNodeValue(node, pos[0], pos[1], pos[2])

    total_faces = 0
    computed_faces = 0
    sti_min = None
    sti_max = None
    sti_sum = 0.0

    for rs in range(rs_count):
        n_faces = reference.GetRsFaceCount(rs)
        out.MakeRs(rs, n_faces, reference.GetRsName(rs), reference.GetRsXmlId(rs))
        print(ui._("   receiver %r: %d faces") % (reference.GetRsName(rs), n_faces))

        for face in range(n_faces):
            total_faces += 1
            if total_faces % PROGRESS_EVERY == 0:
                print("      ... %d / %d faces" % (total_faces, n_faces))

            verts = reference.GetFaceVertices(rs, face)
            out.SetFaceInfo(rs, face, verts[0], verts[1], verts[2], 1)

            # Masking and the threshold need this CELL's absolute band
            # levels, so they are gathered before any MTI is computed.
            cell_factors = None
            if use_masking:
                cell_spl = {}
                for hz in sm.OCTAVE_BANDS_HZ:
                    energy = _face_total_energy(loaded[hz], rs, face)
                    level = (sti_noise.spl_from_energy(energy)
                             if energy is not None else None)
                    if level is not None:
                        cell_spl[hz] = level
                if len(cell_spl) == len(sm.OCTAVE_BANDS_HZ):
                    cell_factors = sm.modulation_reduction_factors(
                        cell_spl, noise_spl_by_band if use_noise else None,
                        True)

            mti_by_band = {}
            usable = True
            for hz in sm.OCTAVE_BANDS_HZ:
                mti = _face_mti(
                    loaded[hz], rs, face, cos_t, sin_t, snr_db_by_band[hz],
                    noise_spl_by_band.get(hz) if use_noise else None,
                    cell_factors[hz] if cell_factors else None)
                if mti is None:
                    usable = False
                    break
                mti_by_band[hz] = mti

            if usable:
                sti = sm.compute_sti(mti_by_band)
                computed_faces += 1
                sti_sum += sti
                sti_min = sti if sti_min is None else min(sti_min, sti)
                sti_max = sti if sti_max is None else max(sti_max, sti)
            else:
                sti = 0.0  # no energy reached this cell

            out.SetFaceEnergy(rs, face, 0, 0, float(sti))

    out_path = _output_path(folder_id)
    if not ls.rsurf_io.Save(out_path, out):
        print(ui._("Failed to write %s") % out_path)
        return

    print(ui._("   computed STI on %d of %d cells") % (computed_faces, total_faces))
    if computed_faces:
        mean = sti_sum / computed_faces
        print("   STI  min %.3f  mean %.3f  max %.3f" % (sti_min, mean, sti_max))
        print("   mean rating: %s (ISO 9921)"
              % sm.iso9921_rating(mean)[0].upper())
    print(ui._("   auditory masking + reception threshold: %s")
          % (u"applied per cell" if use_masking else u"NOT applied"))
    if use_noise:
        print(ui._("   background noise: %s") % noise_description)
        print("      " + "".join("%9d" % hz for hz in sm.OCTAVE_BANDS_HZ) + "   (Hz)")
        print("      " + "".join("%9.1f" % noise_spl_by_band[hz]
                                 for hz in sm.OCTAVE_BANDS_HZ) + "   (noise SPL)")
        print(ui._(
            "      SNR is computed PER CELL from that cell's own speech level, "
            "so it varies across the surface."))
    else:
        print(ui._(
            "   background noise: none - NOISE-FREE map (reverberation only). "
            "This is an upper bound; it will read HIGHER than a punctual "
            "receiver computed with a noise spectrum."))
    print(ui._("   Wrote: %s") % out_path)

    ui.application.sendevent(ui.element(folder_id), ui.idevent.IDEVENT_RELOAD_FOLDER)
    _tag_as_sti_map(folder_id)


def _tag_as_sti_map(folder_id):
    """Marks the freshly-created map element as an STI map so the GUI
    labels it correctly. See the RSTYPE_* notes at the top of this file."""
    target_name = os.path.splitext(OUTPUT_FILENAME)[0]
    try:
        children = ui.element(folder_id).childs()
    except Exception:
        children = []
    for child in children:
        try:
            element = ui.element(child[0])
            if (element.getinfos().get("name") or u"") != target_name:
                continue
            before = _read_rstype(element)
            applied = element.updateentierconfig(RSTYPE_PROPERTY, RSTYPE_STI)
            after = _read_rstype(element)
            print(ui._("   rstype on element %s: before=%s, set=%s -> after=%s")
                  % (child[0], before, "ok" if applied else "FAILED", after))
            if applied and after == RSTYPE_STI:
                print(ui._(
                    "   Tagged as an STI map (rstype=%d). On an I-Simpa build "
                    "that knows the STI record type the legend shows no unit "
                    "and the load menu reads 'Speech Transmision Index (STI)'. "
                    "On older builds (1.3.4 and earlier have no STI support at "
                    "all) the type is unrecognised and the legend falls back to "
                    "'dB' -- IGNORE THE UNIT, the values are dimensionless STI "
                    "in the range 0 to 1.") % RSTYPE_STI)
            else:
                _warn_untagged()
            return
        except Exception:
            continue
    _warn_untagged()


def _read_rstype(element):
    """Reads rstype back so we can see what actually stuck, rather than
    trusting the setter's return value."""
    try:
        return element.getentierconfig(RSTYPE_PROPERTY)
    except Exception as exc:
        return "<unreadable: %s>" % exc


def _warn_untagged():
    print(ui._(
        "   NOTE: could not set rstype=%d on the map element, so the legend "
        "will read seconds (the file itself is written as record type %r "
        "because rsurf has no STI type). The values and colours are still "
        "the STI, 0-1.") % (RSTYPE_STI, OUTPUT_RECORD_TYPE))


def _output_path(folder_id):
    base = ui.e_file(folder_id).buildfullpath()
    folder = base if os.path.isdir(base) else os.path.dirname(base)
    return os.path.join(folder, OUTPUT_FILENAME)


class manager:
    def __init__(self):
        self.compute_map_id = ui.application.register_event(self.OnComputeMap)

    def getmenu(self, typeel, idel, menu):
        try:
            infos = ui.element(idel).getinfos()
        except Exception:
            return False
        name = infos.get("name") or u""
        if name in SURFACE_FOLDER_NAMES:
            menu.insert(0, ())
            menu.insert(0, (ui._("Compute STI Map"), self.compute_map_id))
            return True
        return False

    def OnComputeMap(self, idel):
        try:
            name = ui.element(idel).getinfos().get("name") or u""
        except Exception:
            name = u""
        proceed, noise_spl_by_band, description, use_masking = sti_noise.ask(
            ui._("Compute STI Map for %r") % name)
        if not proceed:
            return
        ComputeSTIMap(idel, noise_spl_by_band, description, use_masking)


ui.application.register_menu_manager(ui.element_type.ELEMENT_TYPE_REPORT_FOLDER, manager())
