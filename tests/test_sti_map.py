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
Tests for the surface-receiver STI map, against a fake libsimpa rsurf that
reproduces the REAL structure captured from a live dump:

  - one .csbin per octave band folder (63 Hz .. 16000 Hz, plus Global)
  - record_type 'SPL_STANDART', 200 steps of 0.01 s, 810 nodes,
    1 surface receiver, 1508 faces
  - per face ~198 SPARSE records, first at step=2, each a
    (time step index, LINEAR energy) pair -- NOT dB

The two properties that would silently corrupt the map if got wrong, and
which these tests pin down:

  1. record index != time step. Records are sparse and start at step 2, so
     times must come from GetFaceTimeStep().
  2. energies are linear, not dB. A stray 10**(dB/10) here would be wrong
     (the punctual-receiver grid IS dB, which is the trap).
"""
import math
import sys
import types

import pytest

# --- fake libsimpa.rsurf, matching the dumped structure -------------------

N_STEPS = 200
TIMESTEP = 0.01
N_NODES = 810
N_FACES = 6           # small stand-in for the real 1508
FIRST_STEP = 2
RECORD_COUNT = 198


def _decay_energy(step, tau):
    t = step * TIMESTEP
    return math.exp(-t / tau)


class FakeRsurfData(object):
    def __init__(self, tau=0.1):
        self.tau = tau
        self._made = None
        self.saved_faces = {}
        self.record_type = "SPL_STANDART"

    # ---- reading ----
    def GetRsCount(self):
        return 1

    def GetNodesCount(self):
        return N_NODES

    def GetTimeStepCount(self):
        return N_STEPS

    def GetTimeStep(self):
        return TIMESTEP

    def getRecordType(self):
        return self.record_type

    def GetRsFaceCount(self, rs):
        return N_FACES

    def GetRsName(self, rs):
        return "Plane Receiver"

    def GetRsXmlId(self, rs):
        return 666

    def GetNodePositionValue(self, i):
        return (float(i), 0.0, 0.0)

    def GetFaceVertices(self, rs, face):
        return (face, face + 1, face + 2)

    def GetFaceRecordCount(self, rs, face):
        return RECORD_COUNT

    def GetFaceTimeStep(self, rs, face, rec):
        # SPARSE: record 0 is step FIRST_STEP, not step 0
        return FIRST_STEP + rec

    def GetFaceEnergy(self, rs, face, rec):
        # LINEAR energy, tiny magnitudes as in the real file
        return 1e-7 * _decay_energy(self.GetFaceTimeStep(rs, face, rec), self.tau)

    def ComputeFaceArea(self, rs, face):
        return 0.1207

    def GetFaceSumEnergy(self, rs, face):
        return sum(self.GetFaceEnergy(rs, face, r) for r in range(RECORD_COUNT))

    # ---- writing ----
    def Make(self, nodes, rs_count, nbtimestep, timestep, record_type):
        self._made = dict(nodes=nodes, rs=rs_count, steps=nbtimestep,
                          dt=timestep, record_type=record_type)

    def MakeRs(self, rs, nfaces, name, xmlid):
        pass

    def SetNodeValue(self, i, x, y, z):
        pass

    def SetFaceInfo(self, rs, face, a, b, c, count):
        pass

    def SetFaceEnergy(self, rs, face, rec, step, energy):
        self.saved_faces[(rs, face)] = energy


SAVED = {}
LOAD_TAUS = {}
RSTYPE_SET = {}


class FakeRsurfIO(object):
    @staticmethod
    def Load(path, data):
        for hz, tau in LOAD_TAUS.items():
            if ("%d Hz" % hz) in path:
                data.tau = tau
                return True
        data.tau = 0.1
        return True

    @staticmethod
    def Save(path, data):
        SAVED[path] = data
        return True


def _install(tmpdir):
    libsimpa = types.ModuleType("libsimpa")
    libsimpa.rsurf_data = FakeRsurfData
    libsimpa.rsurf_io = FakeRsurfIO
    sys.modules["libsimpa"] = libsimpa

    uictrl = types.ModuleType("uictrl")

    class _ET:
        ELEMENT_TYPE_REPORT_FOLDER = "REPFOLDER"

    class _IE:
        IDEVENT_RELOAD_FOLDER = "RELOAD"

    uictrl.element_type = _ET()
    uictrl.idevent = _IE()
    uictrl._ = lambda s: s

    ROOT = 1
    BANDS = [63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]
    band_ids = {hz: 100 + i for i, hz in enumerate(BANDS)}

    MAP_ELEMENT_ID = 900

    class FakeElement(object):
        def __init__(self, eid):
            self.eid = eid

        def getinfos(self):
            if self.eid == ROOT:
                return {"name": "Surface receiver", "label": "Surface receiver",
                        "parentid": -1}
            if self.eid == MAP_ELEMENT_ID:
                return {"name": "STI_map", "label": "STI_map", "parentid": ROOT}
            for hz, bid in band_ids.items():
                if bid == self.eid:
                    return {"name": "%d Hz" % hz, "label": "%d Hz" % hz,
                            "parentid": ROOT}
            raise KeyError(self.eid)

        def childs(self):
            if self.eid == ROOT:
                kids = [(band_ids[hz], "F") for hz in BANDS]
                if SAVED:            # the map element appears after the save
                    kids.append((MAP_ELEMENT_ID, "F"))
                return kids
            return []

        def updateentierconfig(self, name, value):
            RSTYPE_SET[self.eid] = (name, value)
            return True

    class FakeEFile(FakeElement):
        def buildfullpath(self):
            if self.eid == ROOT:
                return str(tmpdir) + "/"
            for hz, bid in band_ids.items():
                if bid == self.eid:
                    import os
                    folder = os.path.join(str(tmpdir), "%d Hz" % hz)
                    if not os.path.isdir(folder):
                        os.makedirs(folder)
                        open(os.path.join(folder, "rs_cut.csbin"), "wb").write(b"x")
                    return folder + "/"
            raise KeyError(self.eid)

    class _App(object):
        pass

    _App.getuserinput = staticmethod(lambda t, m, f: (True, {"Background noise": "None - noise-free STI (reverberation only)"}))
    _App.sendevent = staticmethod(lambda el, ev, p=None: None)
    _App.register_event = staticmethod(lambda fn: "EVT")
    _App.register_menu_manager = staticmethod(lambda t, m: None)

    uictrl.element = FakeElement
    uictrl.e_file = FakeEFile
    uictrl.application = _App()
    sys.modules["uictrl"] = uictrl
    return ROOT


@pytest.fixture
def sti_map_module(tmp_path):
    """Installs a fake uictrl/libsimpa for the map, then RESTORES whatever
    was there before.

    Without the restore this fixture clobbers the shared isimpa_fixtures
    modules for every later test file -- which it did, breaking the
    punctual-receiver tests with a confusing AttributeError. A test file
    must not leak its doubles into the others."""
    saved_modules = {name: sys.modules.get(name)
                     for name in ("uictrl", "libsimpa", "sti_map", "sti_math",
                                  "sti_noise")}
    saved_path = list(sys.path)

    root = _install(tmp_path)
    sys.path.insert(0, _tool_dir('sti_tool'))
    for name in ("sti_map", "sti_math", "sti_noise"):
        sys.modules.pop(name, None)
    import sti_map
    SAVED.clear()
    LOAD_TAUS.clear()
    RSTYPE_SET.clear()
    try:
        yield sti_map, root
    finally:
        sys.path[:] = saved_path
        for name, module in saved_modules.items():
            if module is None:
                sys.modules.pop(name, None)
            else:
                sys.modules[name] = module


def test_sparse_records_are_placed_at_their_own_time_steps(sti_map_module):
    """Record index 0 is time step 2. Using the record index as the time
    would shift the whole echogram earlier and inflate the MTF."""
    sti_map, _root = sti_map_module
    data = FakeRsurfData(tau=0.1)
    assert data.GetFaceTimeStep(0, 0, 0) == FIRST_STEP
    assert data.GetFaceTimeStep(0, 0, 5) == FIRST_STEP + 5

    cos_t, sin_t = sti_map._modulation_tables(N_STEPS, TIMESTEP)
    mti = sti_map._face_mti(data, 0, 0, cos_t, sin_t, 99.0)
    assert mti is not None
    assert 0.0 <= mti <= 1.0


def test_face_mti_matches_the_analytical_mtf_for_an_exponential_decay(sti_map_module):
    """The face echogram here is a pure exponential decay, whose MTF has a
    closed form. This checks the map's per-cell integral independently of
    the punctual-receiver code path."""
    sti_map, _root = sti_map_module
    tau = 0.08
    data = FakeRsurfData(tau=tau)
    cos_t, sin_t = sti_map._modulation_tables(N_STEPS, TIMESTEP)
    mti = sti_map._face_mti(data, 0, 0, cos_t, sin_t, 99.0)

    expected_ti = []
    for f in sm_ref().MODULATION_FREQS_HZ:
        m = 1.0 / math.sqrt(1.0 + (2 * math.pi * f * tau) ** 2)
        expected_ti.append(sm_ref().transmission_index(m))
    expected = sm_ref().mti_for_band(expected_ti)
    assert mti == pytest.approx(expected, abs=0.05)


def sm_ref():
    import sti_math
    return sti_math


def test_linear_energy_is_not_treated_as_db(sti_map_module):
    """Energies are already linear. If a 10**(dB/10) crept in, a tiny value
    like 1e-7 would map to ~1.0 and every face would look identical."""
    sti_map, _root = sti_map_module
    fast = FakeRsurfData(tau=0.02)   # short decay -> high MTI
    slow = FakeRsurfData(tau=1.50)   # long decay  -> low MTI
    cos_t, sin_t = sti_map._modulation_tables(N_STEPS, TIMESTEP)
    mti_fast = sti_map._face_mti(fast, 0, 0, cos_t, sin_t, 99.0)
    mti_slow = sti_map._face_mti(slow, 0, 0, cos_t, sin_t, 99.0)
    assert mti_fast > mti_slow + 0.2


def test_map_writes_one_value_per_face_within_0_and_1(sti_map_module):
    sti_map, root = sti_map_module
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.1
    sti_map.ComputeSTIMap(root)

    assert len(SAVED) == 1
    out = list(SAVED.values())[0]
    assert out._made["record_type"] == sti_map.OUTPUT_RECORD_TYPE
    assert out._made["steps"] == 1
    assert len(out.saved_faces) == N_FACES
    for value in out.saved_faces.values():
        assert 0.0 <= value <= 1.0


def test_more_reverberant_bands_give_a_lower_sti_map(sti_map_module):
    sti_map, root = sti_map_module

    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.03
    sti_map.ComputeSTIMap(root)
    dry = list(list(SAVED.values())[0].saved_faces.values())[0]

    SAVED.clear()
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 1.0
    sti_map.ComputeSTIMap(root)
    wet = list(list(SAVED.values())[0].saved_faces.values())[0]

    assert wet < dry


def test_missing_octave_band_is_reported_and_nothing_written(sti_map_module, capsys):
    sti_map, root = sti_map_module
    original = sti_map._band_folder_map

    def only_low_bands(folder_id):
        full = original(folder_id)
        return {hz: v for hz, v in full.items() if hz <= 1000}

    sti_map._band_folder_map = only_low_bands
    try:
        sti_map.ComputeSTIMap(root)
    finally:
        sti_map._band_folder_map = original

    assert len(SAVED) == 0
    text = capsys.readouterr().out
    assert "Missing octave band" in text


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))


def test_map_element_is_tagged_as_an_sti_map(sti_map_module):
    """rsurf cannot store an STI record type, but the GUI reads the
    'rstype' property off the tree element, not the file. Setting it to 9
    (RECEPTEURS_RECORD_TYPE_STI) is what makes the legend drop the bogus
    'seconds' unit and the load menu read 'Speech Transmision Index'."""
    sti_map, root = sti_map_module
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.1
    sti_map.ComputeSTIMap(root)

    assert sti_map.RSTYPE_STI == 9
    assert RSTYPE_SET, "rstype was never set on the map element"
    name, value = list(RSTYPE_SET.values())[0]
    assert name == "rstype"
    assert value == 9


def test_file_is_still_written_as_a_type_rsurf_accepts(sti_map_module):
    """The FILE must stay one of the four types Make() accepts; only the
    GUI-side property is overridden."""
    sti_map, root = sti_map_module
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.1
    sti_map.ComputeSTIMap(root)
    out = list(SAVED.values())[0]
    assert out._made["record_type"] in ("SPL_STANDART", "SPL_GAIN", "TR", "EDT")
    assert out._made["record_type"] == "TR"


# ---------- background-noise parity with the punctual tool ----------

def test_noise_spectrum_lowers_the_map_versus_noise_free(sti_map_module):
    """The defect this fixes: the map used to be hard-wired noise-free, so
    its MINIMUM STI could come out higher than a punctual receiver
    computed with an NC curve on the same model."""
    sti_map, root = sti_map_module
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.1

    sti_map.ComputeSTIMap(root)
    quiet = list(list(SAVED.values())[0].saved_faces.values())[0]

    SAVED.clear()
    import sti_noise
    noisy_levels = sm_ref().nc_curve_levels(45)
    sti_map.ComputeSTIMap(root, noisy_levels, "NC 45 curve")
    noisy = list(list(SAVED.values())[0].saved_faces.values())[0]

    assert noisy < quiet, (noisy, quiet)


def test_snr_is_derived_per_cell_from_that_cells_own_level(sti_map_module):
    """Speech level varies across the surface, so a single map-wide SNR
    would be wrong everywhere but one spot. A quieter cell must end up
    with a lower STI than a louder one under the same noise."""
    sti_map, _root = sti_map_module
    import sti_noise

    loud = FakeRsurfData(tau=0.1)
    quietened = FakeRsurfData(tau=0.1)
    original = quietened.GetFaceEnergy
    quietened.GetFaceEnergy = lambda rs, f, r: original(rs, f, r) * 1e-3

    cos_t, sin_t = sti_map._modulation_tables(N_STEPS, TIMESTEP)
    noise = 40.0
    mti_loud = sti_map._face_mti(loud, 0, 0, cos_t, sin_t, 99.0, noise)
    mti_quiet = sti_map._face_mti(quietened, 0, 0, cos_t, sin_t, 99.0, noise)
    assert mti_quiet < mti_loud


def test_spl_from_energy_matches_isimpa_display_conversion():
    """I-Simpa shows surface levels as 10*log10(E/1e-12); the per-cell
    speech level must use the same scale as the punctual 'Total' column."""
    import math
    import sti_noise
    assert sti_noise.spl_from_energy(1e-12) == pytest.approx(0.0)
    assert sti_noise.spl_from_energy(1e-6) == pytest.approx(60.0)
    assert sti_noise.spl_from_energy(0.0) is None


# ---------- the map must use the SAME hearing model as the receivers ----------

def test_masking_lowers_the_map_like_it_lowers_a_receiver(sti_map_module):
    """If the map kept the noise-only path it would drift from the punctual
    tool again -- the defect already fixed once for background noise."""
    sti_map, root = sti_map_module
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.1

    sti_map.ComputeSTIMap(root, None, None, False)
    without = list(list(SAVED.values())[0].saved_faces.values())[0]

    SAVED.clear()
    sti_map.ComputeSTIMap(root, None, None, True)
    with_masking = list(list(SAVED.values())[0].saved_faces.values())[0]

    assert with_masking < without, (with_masking, without)


def test_map_reports_whether_the_hearing_model_was_applied(sti_map_module, capsys):
    sti_map, root = sti_map_module
    for hz in sm_ref().OCTAVE_BANDS_HZ:
        LOAD_TAUS[hz] = 0.1
    sti_map.ComputeSTIMap(root, None, None, True)
    assert "applied per cell" in capsys.readouterr().out


def test_noise_dialog_returns_the_masking_choice_too():
    """Both tools read the choice from one shared dialog, so they cannot be
    computed with different assumptions."""
    import sti_noise
    fields = sti_noise.dialog_fields()
    assert sti_noise.LABEL_MASKING in fields
    ok, levels, desc, use_masking = sti_noise.ask("t")
    assert isinstance(use_masking, bool)
