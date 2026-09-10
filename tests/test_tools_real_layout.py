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


"""
Tests both UserScript tools against fixtures replicating the REAL I-Simpa
2.x grid layout (see isimpa_fixtures.py). These replace the earlier
harnesses, which encoded the layout backwards and so passed while the
actual GUI failed.
"""
import math
import sys

import pytest

import isimpa_fixtures as fx

fx.install()

sys.path.insert(0, _userscript_dir())
sys.path.insert(0, _tool_dir('sti_tool'))

import alcons_tool  # noqa: E402
import sti_tool  # noqa: E402
import sti_math as sm  # noqa: E402


def _reset(**kwargs):
    fx.STATE.update({
        "mode": "SPPS",
        "t60": 1.2,
        "nan_bands": set(),
        "omit_bands": set(),
        "echogram_tau": None,
        "user_input": None,
        "receiver_label": "Receiver 1",
    })
    fx.STATE.update(kwargs)
    fx.GABE_WRITES.clear()


STI_INPUT = (True, {"Background noise": "None - noise-free STI (reverberation only)"})
STI_INPUT_NC30 = (True, {"Background noise": "NC 30"})


# ---------------- label / parsing regressions ----------------

def test_rt_column_found_with_real_RT_spelling():
    """The live build labels the column 'RT-30 (s)'; the shipped sample
    script uses 'TR-30 (s)'. Both must be found -- this exact mismatch is
    what made %ALcons fail in the GUI."""
    assert alcons_tool._find_rt_column(fx.AP_HEADER, 30) == fx.AP_RT30_COL
    assert sti_tool._find_rt_column(fx.AP_HEADER, 30) == fx.AP_RT30_COL
    assert alcons_tool._find_rt_column(
        ["", "TR-30 (s)"], 30) == 1
    assert alcons_tool._find_rt_column(["", "Ts (ms)"], 30) == -1


def test_band_label_parsing_excludes_summary_rows():
    for mod in (alcons_tool, sti_tool):
        assert mod._band_label_to_hz("63 Hz") == 63.0
        assert mod._band_label_to_hz("1000 Hz") == 1000.0
        assert mod._band_label_to_hz("16000 Hz") == 16000.0
        assert mod._band_label_to_hz("Global") is None
        assert mod._band_label_to_hz("Average") is None


def test_time_label_parsing_excludes_total_column():
    assert sti_tool._time_label_to_seconds("10.0 ms") == pytest.approx(0.010)
    assert sti_tool._time_label_to_seconds("2000.0 ms") == pytest.approx(2.0)
    assert sti_tool._time_label_to_seconds("Total") is None
    assert sti_tool._time_label_to_seconds("") is None


def test_nan_is_caught_as_invalid():
    for mod in (alcons_tool, sti_tool):
        assert mod._as_finite_float(float("nan")) is None
        assert mod._as_finite_float(float("inf")) is None
        assert mod._as_finite_float("not a number") is None
        assert mod._as_finite_float(1.5) == 1.5


# ---------------- %ALcons ----------------





# ---------------- STI ----------------

def test_sti_spps_reads_real_transposed_echogram():
    """The regression that mattered: with the real layout (bands as rows,
    times as columns) the echogram must actually be found and produce a
    sane STI, where the previous code found nothing."""
    _reset(mode="SPPS", user_input=STI_INPUT)
    echo, missing, decay = sti_tool.GetEchogramPerOctaveBand_SPPS(fx.SOUND_LEVEL_GABE_ID)
    assert echo is not None
    assert missing == []
    for hz in sm.OCTAVE_BANDS_HZ:
        times, h2 = echo[hz]
        assert len(times) == len(fx.TIME_STEPS_MS)
        assert times[0] == pytest.approx(0.010)
        assert times[-1] == pytest.approx(2.0)
        # Leading steps are legitimately ZERO energy: I-Simpa writes -inf dB
        # before the direct sound arrives. Those must be kept (so the time
        # axis stays uniform) but must not be negative, and the band must
        # still carry energy overall.
        assert all(v >= 0 for v in h2)
        assert sum(h2) > 0
        assert h2[0] == 0.0 and h2[fx.STATE["silent_lead_steps"]] > 0

    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    out = fx.written_columns()
    assert 0.0 <= out["STI"] <= 1.0
    for hz in sm.OCTAVE_BANDS_HZ:
        assert 0.0 <= out["MTI %d Hz" % hz] <= 1.0


def test_sti_spps_matches_tcr_for_the_same_decay():
    """Both branches fed the same exponential decay must agree -- the
    cross-check that validates the whole chain end to end."""
    _reset(mode="TCR", user_input=STI_INPUT, t60=1.2)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    sti_tcr = fx.written_columns()["STI"]

    _reset(mode="SPPS", user_input=STI_INPUT, t60=1.2)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    sti_spps = fx.written_columns()["STI"]

    assert sti_spps == pytest.approx(sti_tcr, abs=0.02), (sti_spps, sti_tcr)


def test_sti_tcr_reports_nan_bands_and_refuses(capsys):
    _reset(mode="TCR", user_input=STI_INPUT, nan_bands={4000})
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    assert len(fx.GABE_WRITES) == 0
    text = capsys.readouterr().out
    assert "4000 Hz" in text
    assert "NaN" in text


def test_sti_reports_missing_bands_and_refuses(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT, omit_bands={8000})
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    assert len(fx.GABE_WRITES) == 0
    text = capsys.readouterr().out
    assert "8000 Hz" in text


def test_sti_warns_when_echogram_is_truncated(capsys):
    """A very long decay barely decays within the 2 s window; the tool
    must warn that STI is biased rather than reporting it as-is."""
    _reset(mode="SPPS", user_input=STI_INPUT, echogram_tau=50.0)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "WARNING" in text
    assert "truncated" in text.lower()


def test_sti_does_not_warn_for_a_well_resolved_decay(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT, echogram_tau=0.15)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "WARNING" not in text


# ---------------- textbook / literature provenance ----------------
#
# Cross-checked against Bies & Hansen, "Engineering Noise Control",
# Section 7.11 (Auditorium Design):
#
#   eq (7.129):  %ALcons = 200 * T60^2 * r^2 / V
#                "for the 2000 Hz octave band that is usually used"
#
# This is exactly what alcons_tool implements (and exactly the formula
# specified for this project). NOTE the textbook prints the form WITHOUT
# the directivity factor Q in the denominator; the Peutz/Klein form often
# cited elsewhere is %ALcons = 200*r^2*T60^2/(V*Q). We follow the textbook
# (and the project spec) and use Q only for the critical distance rH.


def test_sti_structure_matches_textbook_description_of_full_sti():
    """Engineering Noise Control 7.11.9.1 describes STI as 14 modulation
    frequencies x 7 octave bands (125 Hz - 8 kHz) = 98 measurements,
    weighted by frequency into a single value -- and RASTI as the reduced
    9-modulation-frequency, 2-band (500 Hz / 2 kHz) shortcut. We implement
    the full 98-point method, not the shortcut."""
    assert len(sm.OCTAVE_BANDS_HZ) == 7
    assert sm.OCTAVE_BANDS_HZ[0] == 125
    assert sm.OCTAVE_BANDS_HZ[-1] == 8000
    assert len(sm.MODULATION_FREQS_HZ) == 14
    assert len(sm.OCTAVE_BANDS_HZ) * len(sm.MODULATION_FREQS_HZ) == 98
    # not the RASTI 2-band shortcut
    assert sm.OCTAVE_BANDS_HZ != [500, 2000]


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))


# ---------- background noise: NC curves and User Defined ----------

def test_noise_mode_dropdown_lists_none_user_defined_then_nc_curves():
    captured = {}

    def fake_getuserinput(title, msg, fields):
        captured["fields"] = dict(fields)
        return (False, {})

    original = fx.uictrl_application_getuserinput_swap(fake_getuserinput)
    try:
        _reset(mode="SPPS")
        sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    finally:
        fx.uictrl_application_getuserinput_swap(original)

    modes = captured["fields"]["Background noise"]
    assert modes[0] == sti_tool.NOISE_MODE_NONE
    assert modes[1] == sti_tool.NOISE_MODE_USER
    assert modes[2:] == ["NC %d" % nc for nc in range(15, 75, 5)]
    assert "NC 15" in modes and "NC 65" in modes and "NC 70" in modes
    assert "RASTI" not in " ".join(modes)


def test_nc_curve_lookup_matches_the_published_table():
    """Cross-checked against two independent published NC tables that
    agreed exactly."""
    nc30 = sm.nc_curve_levels(30)
    assert nc30[125] == 48 and nc30[250] == 41 and nc30[500] == 35
    assert nc30[1000] == 31 and nc30[2000] == 29
    assert nc30[4000] == 28 and nc30[8000] == 27
    nc15 = sm.nc_curve_levels(15)
    assert nc15[125] == 36 and nc15[8000] == 11
    # NC curves slope downward with frequency
    for nc in sm.available_nc_values():
        levels = sm.nc_curve_levels(nc)
        ordered = [levels[hz] for hz in sm.OCTAVE_BANDS_HZ]
        assert all(ordered[i] >= ordered[i + 1] for i in range(len(ordered) - 1))
    # and a higher NC is louder in every band
    for lower, higher in zip(sm.available_nc_values(), sm.available_nc_values()[1:]):
        lo, hi = sm.nc_curve_levels(lower), sm.nc_curve_levels(higher)
        assert all(hi[hz] > lo[hz] for hz in sm.OCTAVE_BANDS_HZ)


def test_nc_noise_lowers_sti_versus_noise_free():
    _reset(mode="SPPS", user_input=STI_INPUT)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    sti_quiet = fx.written_columns()["STI"]

    _reset(mode="SPPS", user_input=STI_INPUT_NC30)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    out = fx.written_columns()
    # the verdict now lives in the ROW LABEL as words, not a 0/1 code
    assert any(k.startswith("Noise included: Yes") for k in out), sorted(out)
    assert any(k.startswith("Background noise: NC 30") for k in out), sorted(out)
    assert out["STI"] <= sti_quiet


def test_noise_free_result_is_flagged_as_such(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "NOISE-FREE" in text
    assert "upper bound" in text


def test_result_reports_both_rating_scales_and_the_source_note(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "ISO 9921" in text
    assert "amplified vs unamplified" in text
    assert "ANSI S3.5" in text


def test_results_file_carries_words_not_numeric_codes():
    """Rating and noise state must read as words in the results
    spreadsheet ('BAD', 'Yes'), not as 0/1 codes."""
    _reset(mode="SPPS", user_input=STI_INPUT)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    keys = list(fx.written_columns().keys())
    joined = " | ".join(keys)

    assert any(k.startswith("Rating (ISO 9921): ") for k in keys), joined
    assert not any("Qualification band" in k for k in keys), joined
    assert "Noise included: No" in joined
    assert "Background noise: none" in joined
    # the old numeric-code labels must be gone
    assert "Rating index" not in joined
    assert "(1=yes)" not in joined


def test_noise_entry_fields_cover_63_hz_to_8k_in_numeric_order():
    """63 Hz is offered even though STI does not use it, and zero-padding
    makes the dialog's alphabetical sort match numeric order."""
    captured = {}

    def fake_getuserinput(title, msg, fields):
        captured["fields"] = dict(fields)
        return (False, {})

    original = fx.uictrl_application_getuserinput_swap(fake_getuserinput)
    try:
        _reset(mode="SPPS")
        sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    finally:
        fx.uictrl_application_getuserinput_swap(original)

    noise_fields = sorted(k for k in captured["fields"] if k.startswith("Noise SPL"))
    freqs = [int(k.split()[2]) for k in noise_fields]
    assert freqs == [63, 125, 250, 500, 1000, 2000, 4000, 8000]
    assert freqs == sorted(freqs), "alphabetical order must match numeric order"



# ---------- "Echogram per source" must never be picked up by accident ----------
#
# With the option ticked, SPPS writes an EXTRA report, also named "Sound
# level", into a sub-folder per source (reportmanager.cpp:620), while the
# COMBINED echogram is still written unconditionally at the receiver
# (reportmanager.cpp:613). I-Simpa's getallelementbytype() recurses
# (element.cpp:306-316), so a first-match search can return ONE source's
# echogram -- which would inflate E_dir badly, because the direct-arrival
# logic marks EVERY source's arrival.

def test_combined_echogram_is_chosen_even_with_per_source_enabled():
    for count in (0, 1, 3):
        _reset(per_source_echograms=count)
        assert sti_tool.FindSoundLevelGabeId(fx.RECEIVER_FOLDER_ID) ==             fx.SOUND_LEVEL_GABE_ID, count
        assert alcons_tool.FindSoundLevelGabeId(fx.RECEIVER_FOLDER_ID) ==             fx.SOUND_LEVEL_GABE_ID, count


def test_a_per_source_report_is_never_returned():
    _reset(per_source_echograms=3)
    got = sti_tool.FindSoundLevelGabeId(fx.RECEIVER_FOLDER_ID)
    assert got not in fx.PER_SOURCE_SOUND_LEVEL_IDS


def test_per_source_echograms_are_counted_so_the_choice_is_visible():
    for count in (0, 1, 3):
        _reset(per_source_echograms=count)
        assert sti_tool.CountPerSourceEchograms(fx.RECEIVER_FOLDER_ID) == count


def test_the_user_is_told_which_echogram_was_used(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT, per_source_echograms=2)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "COMBINED echogram" in text
    assert "NOT used" in text


def test_no_note_when_the_option_is_off(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT, per_source_echograms=0)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    assert "COMBINED echogram" not in capsys.readouterr().out


def test_results_file_labels_and_values_stay_aligned():
    """row_labels and datacol are built in two separate places, so a
    mismatch silently shifts every value against its label."""
    _reset(mode="SPPS", user_input=STI_INPUT)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    _, cols = list(fx.GABE_WRITES.items())[0]
    lengths = {name: len(values) for name, values in cols}
    assert len(set(lengths.values())) == 1, lengths


def test_schroeder_estimate_is_reported_for_an_spps_receiver(capsys):
    _reset(mode="SPPS", user_input=STI_INPUT)
    sti_tool.ComputeSTI(fx.RECEIVER_FOLDER_ID)
    assert "SCHROEDER" in capsys.readouterr().out
    assert any("Schroeder" in k for k in fx.written_columns())
