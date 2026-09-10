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
# Shared background-noise selection for the STI tools.
#
# Both the punctual-receiver STI (sti_tool) and the surface STI map
# (sti_map) must offer the SAME noise choices, or the two disagree on the
# same model -- which is exactly what happened: the map was hard-wired to
# noise-free while the receivers used an NC curve, so the map's MINIMUM
# STI came out higher than the receiver values. This module exists so
# there is one definition rather than two that can drift.
from __future__ import print_function
import os
import sys

import uictrl as ui

_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)
import sti_math as sm

MODE_NONE = u"None - noise-free STI (reverberation only)"
MODE_USER = u"User Defined"

DEFAULT_NOISE_SPL_DB = 30.0
# Expressed as a very large SNR rather than a special case, so the same
# code path runs either way: 10^(-99/10) ~ 1e-10 makes the noise term in
# m'(F) = m(F)/(1+10^(-SNR/10)) vanish.
NOISE_FREE_SNR_DB = 99.0

# 63 Hz is offered even though STI uses 125 Hz - 8 kHz, so the spectrum
# matches the familiar 8-band NC table. Zero-padding makes the dialog's
# alphabetical field order match numeric order.
ENTRY_BANDS_HZ = [63] + list(sm.OCTAVE_BANDS_HZ)

LABEL_MODE = u"Background noise"

# The IEC hearing model is offered wherever noise is, so the map and the
# punctual receivers can never be computed with different assumptions.
LABEL_MASKING = u"Auditory masking + reception threshold"
MASKING_ON = u"On - full IEC hearing model (recommended)"
MASKING_OFF = u"Off - reverberation and ambient noise only"
MASKING_MODES = [MASKING_ON, MASKING_OFF]


def nc_mode_label(nc):
    return u"NC %d" % nc


MODES = ([MODE_NONE, MODE_USER]
         + [nc_mode_label(nc) for nc in sm.available_nc_values()])


def noise_label(hz):
    return u"Noise SPL %04d Hz (dB)" % hz


def dialog_fields():
    fields = {LABEL_MODE: list(MODES),
              LABEL_MASKING: list(MASKING_MODES)}
    for hz in ENTRY_BANDS_HZ:
        fields[noise_label(hz)] = "%g" % DEFAULT_NOISE_SPL_DB
    return fields


DIALOG_HELP = (
    "Choose the background noise: none, an NC curve, or User Defined.\n\n"
    "IMPORTANT: the 'Noise SPL' fields below are used ONLY when "
    "'User Defined' is selected. Picking an NC curve does NOT rewrite them "
    "(this dialog cannot update its own fields) -- the NC levels are applied "
    "internally and the exact spectrum used is reported with the result.\n\n"
    "The per-band SNR is computed from the speech level predicted at the "
    "receiver minus the background noise level."
)


def resolve(answers):
    """Turns the dialog answers into (ok, noise_spl_by_band, description).

    noise_spl_by_band is empty for the noise-free case, which callers can
    test with `if noise_spl_by_band:`."""
    mode = answers.get(LABEL_MODE, MODE_NONE)

    if mode == MODE_NONE:
        return True, {}, mode

    if mode == MODE_USER:
        levels = {}
        try:
            for hz in sm.OCTAVE_BANDS_HZ:
                levels[hz] = float(answers[noise_label(hz)])
        except (ValueError, KeyError):
            print(ui._("Invalid background noise level."))
            return False, {}, mode
        return True, levels, u"user defined spectrum"

    for nc in sm.available_nc_values():
        if mode == nc_mode_label(nc):
            return True, sm.nc_curve_levels(nc), u"NC %d curve" % nc

    print(ui._("Unrecognised background noise selection: %s") % mode)
    return False, {}, mode


def ask(title):
    """Shows the dialog. Returns
    (proceed, noise_spl_by_band, description, use_masking).
    proceed is False if the user cancelled or gave bad input."""
    result = ui.application.getuserinput(title, ui._(DIALOG_HELP), dialog_fields())
    if not result or not result[0]:
        return False, {}, u"cancelled", True
    ok, levels, description = resolve(result[1])
    use_masking = result[1].get(LABEL_MASKING, MASKING_ON) == MASKING_ON
    return ok, levels, description, use_masking


def spl_from_energy(total_energy):
    """Surface-receiver face energy -> SPL in dB.

    I-Simpa itself displays surface-receiver levels as
    10*log10(E / 1e-12)  (Recepteurs_surfacique.cpp, where the auto-scale
    converts minfacetread/maxfacetread for SPL maps), so the same
    conversion is used here to put the per-cell speech level on the same
    scale as the punctual receiver's 'Total' column."""
    if total_energy <= 0:
        return None
    import math
    return 10.0 * math.log10(total_energy / 1e-12)
