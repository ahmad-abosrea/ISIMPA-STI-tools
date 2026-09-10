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
"""
Tests for alcons_tool, which now implements ONE method: the energetic form

    %ALcons = 8.9 * T60 * E_rev / (13.82 * E_dir)

(Bies & Hansen, "Engineering Noise Control", eq 7.130, after Bistafa &
Bradley 2000). The classical single-source form (eq 7.129) was removed
deliberately: it needs an r and a Q, which for a multi-source project is
both ambiguous and a way for a user to silently contradict the directivity
actually configured on the sources.

The classical form is still used HERE as an independent reference to
verify the energetic constant, even though the tool no longer offers it.
"""
import math
import sys

import pytest

import isimpa_fixtures as fx

fx.install()
import alcons_tool  # noqa: E402


def _classical_reference(t60, r, volume, q):
    """eq (7.129) with Q -- reference only; not part of the tool."""
    return 200.0 * (r ** 2) * (t60 ** 2) / (volume * q)


def _rev_over_dir(r, t60, volume, q):
    """Steady-state: E_dir = W*Q/(4*pi*r^2*c), E_rev = 4W/(R*c)
       => E_rev/E_dir = 16*pi*r^2/(Q*R), Sabine R ~= A = 0.161V/T60."""
    absorption = 0.161 * volume / t60
    return 16.0 * math.pi * r * r / (q * absorption)


# ---------- is the constant right? three independent checks ----------

@pytest.mark.parametrize("r,t60,volume,q", [
    (1.0, 1.0, 500.0, 1.0),
    (2.0, 1.0, 500.0, 5.0),
    (2.0, 2.0, 2000.0, 2.0),
    (5.0, 3.5, 5000.0, 10.0),
    (3.0, 1.5, 1200.0, 4.0),
])
def test_energetic_form_reproduces_the_classical_peutz_coefficient(r, t60, volume, q):
    """Substituting the steady-state direct/reverberant ratio into
    eq (7.130) must reproduce the classical 200 coefficient (it gives
    201.1, a 0.5% match). This is what justifies the 8.9/13.82 constant."""
    ratio = _rev_over_dir(r, t60, volume, q)
    energetic = alcons_tool.AlconsRatioTerm(t60, e_dir=1.0, e_rev=ratio)
    assert energetic == pytest.approx(_classical_reference(t60, r, volume, q), rel=0.01)


@pytest.mark.parametrize("t60,volume,q", [
    (1.0, 500.0, 1.0), (2.0, 2000.0, 2.0), (3.5, 5000.0, 10.0), (1.2, 200.0, 1.0),
])
def test_agrees_with_classical_at_the_critical_distance(t60, volume, q):
    """At the critical distance rH, E_dir = E_rev BY DEFINITION, so both
    formulations must agree there. Classical at rH gives 0.65*T60; the
    energetic form with ratio = 1 gives 0.644*T60."""
    rh = 0.057 * math.sqrt(q * volume / t60)
    classical_at_rh = _classical_reference(t60, rh, volume, q)
    energetic_at_unity_ratio = alcons_tool.AlconsRatioTerm(t60, e_dir=1.0, e_rev=1.0)
    assert energetic_at_unity_ratio == pytest.approx(classical_at_rh, rel=0.02)


def test_the_9x_variant_would_be_wrong_by_ln_10_to_the_6():
    """Guards against 'simplifying' to 9*(E_rev/E_dir)*T60, which drops
    the 13.82 = ln(10^6) divisor and overestimates ~13.8x. If someone
    deliberately sets the constants to the 9x scaling this test is the
    thing that should make them stop and think."""
    assert alcons_tool.ALCONS_NUMERATOR == 8.9
    assert alcons_tool.ALCONS_DENOMINATOR == pytest.approx(13.82)
    assert alcons_tool.ALCONS_DENOMINATOR == pytest.approx(math.log(10 ** 6), rel=0.001)

    t60, ratio = 2.0, 3.0
    correct = alcons_tool.AlconsRatioTerm(t60, e_dir=1.0, e_rev=ratio)
    naive = 9.0 * ratio * t60
    assert naive / correct == pytest.approx(13.82 * 9.0 / 8.9, rel=0.01)


def test_formula_text_reports_the_active_constants():
    """The active formula is printed with every result so the scaling in
    use can never be ambiguous."""
    text = alcons_tool.formula_text()
    assert "8.9" in text and "13.82" in text


def _reset(**kwargs):
    fx.STATE.update({
        "mode": "SPPS", "t60": 1.2, "nan_bands": set(), "omit_bands": set(),
        "echogram_tau": None, "receiver_label": "Receiver 1",
        "temperature_c": 20.0, "has_environment": True,
        "user_input": (True, {"Algorithm": alcons_tool.METHOD_BISTAFA}),
    })
    fx.STATE.update(kwargs)
    fx.GABE_WRITES.clear()


# ---------- the tool exposes no Q / V / r / source knobs any more ----------

def test_classical_entry_points_are_gone():
    for gone in ("AlconsClassical", "METHOD_CLASSICAL", "METHOD_ENERGETIC"):
        assert not hasattr(alcons_tool, gone), "%s should have been removed" % gone


def sti_noise_fields():
    import sti_noise
    return sti_noise.dialog_fields().keys()


def test_dialog_only_offers_an_algorithm_choice():
    """The only thing asked of the user is WHICH PUBLISHED ALGORITHM to
    report -- no physical quantity, so nothing can contradict the model."""
    captured = {}

    def fake_getuserinput(title, msg, fields):
        captured["fields"] = dict(fields)
        return (False, {})

    original = fx.uictrl_application_getuserinput_swap(fake_getuserinput)
    try:
        _reset()
        alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    finally:
        fx.uictrl_application_getuserinput_swap(original)

    keys = set(captured["fields"].keys())
    assert captured["fields"]["Algorithm"] == [
        alcons_tool.METHOD_BISTAFA, alcons_tool.METHOD_PEUTZ_LIMIT]

    # Background noise IS asked for -- it is genuinely external to the model
    # (the room file says nothing about HVAC noise), and it feeds only the
    # STI-derived cross-check. What must stay absent is any quantity the
    # model already knows, which a user could set to contradict it.
    for forbidden in ("Q", "Directivity", "Volume", "Distance", "r", "Source"):
        assert not any(forbidden.lower() == k.lower() for k in keys), keys
    assert keys - {"Algorithm"} <= set(sti_noise_fields()), keys


# ---------- the two published constants are one formula ----------

def test_the_9_is_the_ceiling_not_the_ratio_multiplier():
    """Peutz: "without direct sound %ALcons is limited to 9*T60". The ratio
    term reaches that ceiling at E_rev/E_dir = 13.85 (r = 3.72*Dc, i.e.
    D/R = -11.4 dB, matching the published "D/R = -11 dB at 3.5*Dc"). So a
    ratio form C*(ratio)*T60 needs C = 9/13.85 = 0.650 -- which is exactly
    Bistafa & Bradley's 8.9/13.82 = 0.644."""
    assert alcons_tool.ALCONS_CEILING_FACTOR == 9.0
    c_used = alcons_tool.ALCONS_NUMERATOR / alcons_tool.ALCONS_DENOMINATOR
    assert c_used == pytest.approx(9.0 / 13.85, rel=0.02)

    t60 = 2.0
    ceiling = alcons_tool.AlconsPeutzCeiling(t60)
    assert alcons_tool.AlconsRatioTerm(t60, 1.0, 13.85) == pytest.approx(ceiling, rel=0.02)
    assert alcons_tool.AlconsRatioTerm(t60, 1.0, 1.0) < 0.1 * ceiling


def test_result_can_never_exceed_the_peutz_ceiling():
    """The saturation cap is what makes impossible (>100%) values
    unreachable, however reverberant the room."""
    t60 = 3.8
    ceiling = alcons_tool.AlconsPeutzCeiling(t60)
    for ratio in [1, 5, 13.85, 50, 500, 100000]:
        value, saturated = alcons_tool.AlconsEnergetic(t60, 1.0, ratio)
        assert value <= ceiling + 1e-9
        if ratio > 20:
            assert saturated is True
    assert ceiling < 100.0


def test_zero_absorption_room_lands_exactly_on_the_peutz_limit():
    """A room with no absorption has essentially no direct-to-reverberant
    ratio left, which is precisely Peutz's "without direct sound" case, so
    the answer must be exactly 9*T60."""
    t60 = 3.8
    value, saturated = alcons_tool.AlconsEnergetic(t60, e_dir=1.0, e_rev=1e6)
    assert saturated is True
    assert value == pytest.approx(9.0 * t60)


# ---------- eq (7.131) direct-to-reverberant rating ----------

def test_direct_to_reverberant_rating_scale():
    """> -3 dB excellent ... < -15 dB very poor, in 3 dB steps."""
    assert alcons_tool.rate_direct_to_reverberant(0.0) == "excellent"
    assert alcons_tool.rate_direct_to_reverberant(-2.9) == "excellent"
    assert alcons_tool.rate_direct_to_reverberant(-4.0) == "very good"
    assert alcons_tool.rate_direct_to_reverberant(-7.0) == "good"
    assert alcons_tool.rate_direct_to_reverberant(-10.0) == "fair"
    assert alcons_tool.rate_direct_to_reverberant(-13.0) == "poor"
    assert alcons_tool.rate_direct_to_reverberant(-16.0) == "very poor"
    assert alcons_tool.rate_direct_to_reverberant(-40.0) == "very poor"


# ---------- speed of sound from the project temperature ----------

def test_celerity_matches_isimpa_own_formula():
    """I-Simpa: c_son(K) = 343.2 * sqrt(K / 293.15), K = degC + 273.15
    (Celerite_du_son.cpp + calculsPropagation.h Kref). Matching it exactly
    keeps our arrival times consistent with how the engine propagated the
    particles."""
    assert alcons_tool.speed_of_sound_from_temperature(20.0) == pytest.approx(343.2)
    for t_c in (0.0, 10.0, 20.0, 25.0, 35.0):
        expected = 343.2 * math.sqrt((t_c + 273.15) / 293.15)
        assert alcons_tool.speed_of_sound_from_temperature(t_c) == pytest.approx(expected)
    # warmer air is faster
    assert (alcons_tool.speed_of_sound_from_temperature(30.0)
            > alcons_tool.speed_of_sound_from_temperature(10.0))


def test_temperature_is_read_from_the_project_environment():
    _reset(temperature_c=30.0)
    temp, found = alcons_tool.GetProjectTemperature()
    assert found is True
    assert temp == pytest.approx(30.0)


def test_falls_back_and_warns_when_environment_is_missing(capsys):
    _reset(has_environment=False)
    temp, found = alcons_tool.GetProjectTemperature()
    assert found is False
    assert temp == pytest.approx(20.0)

    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "WARNING" in text and "Environment" in text


def test_reported_celerity_follows_project_temperature(capsys):
    _reset(temperature_c=30.0)
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    expected_c = alcons_tool.speed_of_sound_from_temperature(30.0)
    assert "%.1f m/s" % expected_c in text
    assert "30.0 degC" in text


# ---------- the direct/reverberant split ----------

def test_split_assigns_each_source_arrival_to_its_own_step():
    times = [0.010 * i for i in range(1, 21)]
    energies = [1.0] * len(times)
    c = 343.8
    # 5.0 m -> 14.5 ms lands in step (10,20]; 20.0 m -> 58.2 ms in (50,60]
    e_dir, e_rev, direct_times, dt = alcons_tool.SplitDirectReverberant(
        times, energies, [5.0 / c, 20.0 / c])
    assert dt == pytest.approx(0.010)
    assert direct_times == pytest.approx([0.020, 0.060])
    assert e_dir == pytest.approx(2.0)
    assert e_rev == pytest.approx(len(times) - 2.0)


def test_split_counts_a_shared_step_once():
    times = [0.010 * i for i in range(1, 11)]
    energies = [1.0] * len(times)
    c = 343.8
    _e_dir, _e_rev, direct_times, _ = alcons_tool.SplitDirectReverberant(
        times, energies, [6.9 / c, 7.0 / c])
    assert len(direct_times) == 1


def test_split_energy_is_conserved():
    times = [0.010 * i for i in range(1, 31)]
    energies = [float(i) for i in range(1, 31)]
    e_dir, e_rev, _, _ = alcons_tool.SplitDirectReverberant(times, energies, [0.05, 0.12])
    assert e_dir + e_rev == pytest.approx(sum(energies))


def test_more_sources_raise_direct_energy_and_lower_alcons():
    """Adding a second source contributes its own direct arrival, raising
    E_dir and improving (lowering) %ALcons -- the multi-source behaviour
    the energetic method exists for."""
    times = [0.010 * i for i in range(1, 21)]
    energies = [1.0] * len(times)
    c = 343.8
    one = alcons_tool.SplitDirectReverberant(times, energies, [5.0 / c])
    two = alcons_tool.SplitDirectReverberant(times, energies, [5.0 / c, 20.0 / c])
    a_one = alcons_tool.AlconsEnergetic(1.2, one[0], one[1])
    a_two = alcons_tool.AlconsEnergetic(1.2, two[0], two[1])
    assert two[0] > one[0]
    assert a_two < a_one


# ---------- end to end through the fake GUI ----------

def _reset(**kwargs):
    fx.STATE.update({
        "mode": "SPPS", "t60": 1.2, "nan_bands": set(), "omit_bands": set(),
        "echogram_tau": None, "receiver_label": "Receiver 1",
        "temperature_c": 20.0, "has_environment": True,
        "user_input": (True, {"Algorithm": alcons_tool.METHOD_BISTAFA}),
    })
    fx.STATE.update(kwargs)
    fx.GABE_WRITES.clear()


def test_runs_end_to_end_and_is_self_consistent():
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    out = fx.written_columns()
    assert out["T60 2kHz (s)"] == pytest.approx(1.2)
    assert out["Number of ACTIVE sources"] == 1.0
    assert out["E_rev / E_dir"] > 0
    expected, _sat = alcons_tool.AlconsEnergetic(1.2, 1.0, out["E_rev / E_dir"])
    assert out["%ALcons energy-based (%) - diagnostic only"] == pytest.approx(expected, rel=1e-6)
    # no Q / volume / r columns should survive
    for gone in ("Q used", "r source-receiver (m)", "rH critical distance (m)"):
        assert gone not in out


def test_refuses_and_explains_when_2khz_rt_is_nan(capsys):
    _reset(nan_bands={2000})
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    assert len(fx.GABE_WRITES) == 0
    text = capsys.readouterr().out
    assert "NaN" in text and "particle" in text.lower()


def test_warns_when_time_step_is_too_coarse(capsys):
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "WARNING" in text and "time step" in text.lower()


def test_prints_every_source_and_its_arrival(capsys):
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    assert "SpeakerA" in text
    assert "direct arrival" in text
    assert "directivities as configured" in text


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))


def test_peutz_limit_algorithm_reports_the_ceiling():
    _reset(user_input=(True, {"Algorithm": alcons_tool.METHOD_PEUTZ_LIMIT}))
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    out = fx.written_columns()
    assert out["%ALcons energy-based (%) - diagnostic only"] == pytest.approx(9.0 * 1.2)
    # the "Peutz limit" row was removed as redundant (it is just 9 * T60),
    # so the ceiling is checked against the function that produces it
    assert alcons_tool.AlconsPeutzCeiling(1.2) == pytest.approx(9.0 * 1.2)


def test_both_algorithms_and_dr_bands_are_reported(capsys):
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    out = fx.written_columns()
    for key in ("STI (same receiver)",
                "%ALcons from STI, Becker (%) - REPORT THIS",
                "%ALcons energy-based (%) - diagnostic only",
                "Agreement ratio (Becker / energy)",
                "Saturated at Peutz limit (1=yes)",
                "D/R 500 Hz (dB)", "D/R 1000 Hz (dB)", "D/R 2000 Hz (dB)"):
        assert key in out
    text = capsys.readouterr().out
    assert "eq (7.131)" in text
    assert "Algorithm:" in text


# ---------- regressions from live GUI failures ----------

def test_output_path_uses_e_file_not_element():
    """buildfullpath() exists only on uictrl.e_file. The fakes used to put
    it on element too, so `receiver.buildfullpath()` passed in tests and
    raised AttributeError in the real GUI. element must NOT have it."""
    assert not hasattr(fx.uictrl_element_class()(fx.RECEIVER_FOLDER_ID), "buildfullpath")
    path = alcons_tool.receiver_output_path(fx.RECEIVER_FOLDER_ID, "ALcons.gabe")
    assert path.endswith("ALcons.gabe")


def test_minus_inf_steps_are_zero_energy_not_dropped():
    """I-Simpa writes -inf dB for steps before the direct sound arrives.
    That is zero energy, not missing data: the step must be kept so the
    time axis stays uniform."""
    assert alcons_tool._db_to_linear_energy(float("-inf")) == 0.0
    assert alcons_tool._db_to_linear_energy(0.0) == pytest.approx(1.0)
    assert alcons_tool._db_to_linear_energy(float("nan")) is None
    assert alcons_tool._db_to_linear_energy("x") is None

    _reset()
    times, energies, _hz = alcons_tool.GetEchogramForBand(
        fx.SOUND_LEVEL_GABE_ID, alcons_tool.TARGET_BAND_HZ)
    lead = fx.STATE["silent_lead_steps"]
    assert len(times) == len(fx.TIME_STEPS_MS)      # nothing dropped
    assert energies[0] == 0.0
    assert energies[lead] > 0.0


def test_direct_step_skips_leading_silence():
    """If the geometric arrival lands in a step that still holds zero
    energy, the direct sound is the first step that actually carries
    energy -- otherwise E_dir would be 0 and the result undefined."""
    times = [0.010 * i for i in range(1, 11)]
    energies = [0.0, 0.0, 5.0, 4.0, 3.0, 2.0, 1.0, 1.0, 1.0, 1.0]
    # arrival computed at 12 ms -> step (10,20] which is still silent
    e_dir, e_rev, direct_times, _dt = alcons_tool.SplitDirectReverberant(
        times, energies, [0.012])
    assert direct_times == pytest.approx([0.030])
    assert e_dir == pytest.approx(5.0)
    assert e_rev == pytest.approx(sum(energies) - 5.0)


def test_runs_end_to_end_with_leading_silence():
    _reset(silent_lead_steps=3)
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    out = fx.written_columns()
    assert out["%ALcons energy-based (%) - diagnostic only"] > 0
    assert out["E_rev / E_dir"] > 0


# ---------- the STI-derived cross-check (secondary, empirical) ----------

def test_becker_conversion_reproduces_the_published_correspondence_table():
    """The conversion is empirical and weakly sourced, so what is pinned is
    that it reproduces the PUBLISHED STI/%ALcons table it claims to."""
    import sti_math as sm
    for sti, published in ((0.30, 33.0), (0.45, 15.0), (0.60, 7.0), (0.75, 3.0)):
        got, _trusted = sm.alcons_from_sti(sti)
        assert abs(got - published) < 0.7, (sti, got, published)


def test_becker_round_trips_with_its_own_inverse():
    import sti_math as sm
    for sti in (0.25, 0.40, 0.55, 0.70, 0.85):
        back = sm.sti_from_alcons(sm.alcons_from_sti(sti)[0])
        assert abs(back - sti) < 0.001, (sti, back)


def test_becker_is_capped_and_flagged_below_its_valid_range():
    """It runs to 170 % at STI 0, where %ALcons is conventionally capped at
    100 %. Below STI 0.20 it is extrapolation and must say so."""
    import sti_math as sm
    value, trusted = sm.alcons_from_sti(0.0)
    assert value == 100.0 and trusted is False
    assert sm.alcons_from_sti(0.19)[1] is False
    assert sm.alcons_from_sti(0.21)[1] is True


def test_energy_result_is_still_the_primary_and_is_noise_independent():
    """The energy formula must not start depending on background noise --
    that is precisely the limitation the cross-check exists to expose."""
    quiet = alcons_tool.AlconsEnergetic(2.0, 1.0, 5.0)
    assert quiet == alcons_tool.AlconsEnergetic(2.0, 1.0, 5.0)
    import inspect
    src = inspect.getsource(alcons_tool.AlconsEnergetic)
    assert "noise" not in src.lower()


def test_divergence_threshold_is_a_factor_of_two():
    assert alcons_tool.DIVERGENCE_FACTOR == 2.0


# ---------- only ACTIVE sources may contribute a direct arrival ----------
#
# A source with "Active source" unticked is not written into the
# calculation XML at all (e_scene_sources_source.h:104-111), so SPPS never
# traces it and it radiates nothing. Marking a direct-sound step at r/c for
# such a source moves a step that holds only OTHER sources' reflections out
# of E_rev and into E_dir, making %ALcons look far better than it is -- and
# the error grows with every disabled source.

def test_inactive_sources_are_excluded_from_the_direct_split():
    fx.STATE["inactive_sources"] = ("SpeakerA",)
    try:
        assert alcons_tool.GetSceneSourcePositions() == {}
        assert alcons_tool.GetInactiveSceneSourceNames() == ["SpeakerA"]
    finally:
        fx.STATE.pop("inactive_sources", None)


def test_active_sources_are_kept():
    fx.STATE.pop("inactive_sources", None)
    assert "SpeakerA" in alcons_tool.GetSceneSourcePositions()
    assert alcons_tool.GetInactiveSceneSourceNames() == []


def test_the_enable_flag_is_read_from_the_properties_child_not_the_source():
    """Element::GetBoolConfig only sees DIRECT children, so the flag must
    be found by walking the source's children. Asking the source element
    itself must not be what decides it."""
    fx.STATE["inactive_sources"] = ("SpeakerA",)
    try:
        assert alcons_tool.SourceIsActive(1) is False
    finally:
        fx.STATE.pop("inactive_sources", None)
    assert alcons_tool.SourceIsActive(1) is True


def test_a_source_with_no_flag_counts_as_active():
    """Projects whose sources carry no 'enable' property must behave as
    before rather than silently losing every source."""
    assert alcons_tool.SourceIsActive(999999) is True


def test_all_sources_disabled_is_reported_clearly(capsys):
    fx.STATE["inactive_sources"] = ("SpeakerA",)
    try:
        _reset()
        alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
        text = capsys.readouterr().out
        assert "Active source" in text and "Activate at least one" in text
    finally:
        fx.STATE.pop("inactive_sources", None)


# ---------- the trimmed results table ----------

def _labels_after_run():
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    return list(fx.written_columns().keys())


def test_redundant_rows_are_gone():
    """Three rows carried no information of their own: a duplicate of the
    D/R 2 kHz row, 9*T60 restated, and the uncapped ratio term which only
    differs from %ALcons when the saturation flag is already set."""
    labels = _labels_after_run()
    joined = " | ".join(labels)
    assert "Direct-to-reverberant 2kHz" not in joined, joined
    assert "Peutz limit 9*T60 (%)" not in joined, joined
    assert "ratio term, uncapped" not in joined, joined


def test_the_sti_derived_value_is_listed_first_and_marked_as_the_one_to_report():
    labels = _labels_after_run()
    assert labels[0].startswith("STI"), labels[0]
    assert "REPORT THIS" in labels[1], labels[1]
    assert "Becker" in labels[1]


def test_the_energy_value_is_kept_but_marked_diagnostic():
    """It must stay visible -- it is the only independent estimate -- but
    it must not read as the headline for a multi-source model."""
    labels = _labels_after_run()
    energy = [l for l in labels if "energy-based" in l]
    assert len(energy) == 1, labels
    assert "diagnostic only" in energy[0]


def test_the_inputs_that_drive_the_known_failure_modes_are_retained():
    """Source count and time step stay in the file: together they are what
    makes the energy-based split unreliable, so a saved result must carry
    enough to judge it later."""
    joined = " | ".join(_labels_after_run())
    for needed in ("Number of ACTIVE sources", "Echogram time step (ms)",
                   "E_rev / E_dir", "T60 2kHz (s)",
                   "Agreement ratio (Becker / energy)"):
        assert needed in joined, (needed, joined)


def test_per_band_direct_to_reverberant_rows_survive():
    joined = " | ".join(_labels_after_run())
    for hz in (500, 1000, 2000):
        assert "D/R %d Hz (dB)" % hz in joined, joined


def test_labels_and_values_stay_aligned():
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    _, cols = list(fx.GABE_WRITES.items())[0]
    lengths = {name: len(values) for name, values in cols}
    assert len(set(lengths.values())) == 1, lengths


def test_console_guidance_agrees_with_the_results_table(capsys):
    """The report and the results file must not name different winners.
    They briefly did: the table said "REPORT THIS" on the STI-derived row
    while the console said the energy-based one was primary."""
    _reset()
    alcons_tool.ComputeALcons(fx.RECEIVER_FOLDER_ID)
    text = capsys.readouterr().out
    labels = " | ".join(fx.written_columns().keys())

    assert "REPORT THIS" in labels
    reported = [l for l in fx.written_columns() if "REPORT THIS" in l][0]
    assert "from STI" in reported                       # the table's choice
    assert "WHICH TO REPORT" in text                    # console explains it
    assert "report this one" in text                    # and agrees
    assert "the primary result above is the energy-based one" not in text
