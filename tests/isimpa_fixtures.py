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
Shared fake `uictrl` / `libsimpa` built to match the REAL I-Simpa 2.x data
layout, captured from a live debug dump of the user's project (see
debug_tool). The previous version of these harnesses encoded the layout
BACKWARDS (time steps as rows, bands as columns) and therefore passed
while the real scripts failed in the GUI -- these fixtures exist so that
class of error cannot silently repeat.

Ground truth from the live dump:

  'Sound level' report (202 columns x 11 rows):
    row[0] = ['', '10.0 ms', '20.0 ms', ..., '2000.0 ms', 'Total']
    row[i] = ['63 Hz', <200 dB values>, <total dB>]
    -> ROWS = frequency bands, COLUMNS = time steps.

  'Acoustic parameters' report (10 columns x 12 rows):
    row[0] = ['', 'Sound level (dB)', 'Sound level (dBA)', 'C-50 (dB)',
              'C-80 (dB)', 'D-50 (%)', 'Ts (ms)', 'RT-30 (s)', 'EDT (s)',
              'ST (dB)']
    row[i] = ['63 Hz', <9 values>]
    last rows = 'Global' / 'Average' summaries.
    -> note 'RT-30 (s)', NOT 'TR-30 (s)'.
"""
import math
import sys
import types

# Bands actually present in the user's project (63 Hz .. 16 kHz).
PROJECT_BANDS_HZ = [63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000]

# Echogram time grid from the dump: 10 ms steps out to 2000 ms.
TIME_STEPS_MS = [10.0 * i for i in range(1, 201)]

AP_HEADER = [
    "", "Sound level (dB)", "Sound level (dBA)", "C-50 (dB)", "C-80 (dB)",
    "D-50 (%)", "Ts (ms)", "RT-30 (s)", "EDT (s)", "ST (dB)",
]
AP_RT30_COL = AP_HEADER.index("RT-30 (s)")

PUNCTUAL_RECEIVERS_GROUP_ID = 6270
MODE_ANCESTOR_ID = 6260
RECEIVER_FOLDER_ID = 6272
SOUND_LEVEL_GABE_ID = 6276

# "Echogram per source" makes SPPS write an EXTRA report, also named
# "Sound level", into a sub-folder per source. These model those. The
# tools must always read SOUND_LEVEL_GABE_ID (the combined echogram) and
# never one of these.
PER_SOURCE_FOLDER_IDS = [7001, 7002, 7003]
PER_SOURCE_SOUND_LEVEL_IDS = [7101, 7102, 7103]


def _per_source_count():
    return int(STATE.get("per_source_echograms", 0))
ACOUSTIC_PARAMS_ID = 6420
ENVIRONMENT_ID = 6100

# --- knobs the tests set to shape the fixture data ---
STATE = {
    "mode": "SPPS",
    "t60": 1.2,                  # uniform T60 used for both branches
    "nan_bands": set(),          # bands whose RT-30 comes back NaN
    "omit_bands": set(),         # bands entirely absent from the project
    "echogram_tau": None,        # if set, echogram decays with this tau
    "silent_lead_steps": 2,      # steps of -inf dB before the direct sound
    "user_input": None,
    "receiver_label": "Receiver 1",
    "temperature_c": 20.0,      # project Environment temperature
    "has_environment": True,    # set False to test the fallback path
}
GABE_WRITES = {}
PRINTED = []


def _bands():
    return [hz for hz in PROJECT_BANDS_HZ if hz not in STATE["omit_bands"]]


def _make_sound_level_grid():
    """ROWS = bands, COLUMNS = time steps (the real layout)."""
    header = [""] + ["%.1f ms" % ms for ms in TIME_STEPS_MS] + ["Total"]
    rows = [header]
    tau = STATE["echogram_tau"]
    if tau is None:
        tau = STATE["t60"] / 13.8
    lead = STATE.get("silent_lead_steps", 0)
    for hz in _bands():
        values = []
        for i, ms in enumerate(TIME_STEPS_MS):
            if i < lead:
                # I-Simpa writes -inf dB for steps with no energy yet
                values.append(float("-inf"))
                continue
            t = ms / 1000.0
            h2 = math.exp(-t / tau)
            values.append(10 * math.log10(h2))
        finite = [v for v in values if v != float("-inf")]
        total = 10 * math.log10(sum(10 ** (v / 10.0) for v in finite))
        rows.append(["%d Hz" % hz] + values + [total])
    rows.append(["Global"] + [0.0] * (len(TIME_STEPS_MS) + 1))
    return rows


def _make_acoustic_params_grid():
    """ROWS = bands, COLUMNS = parameters (the real layout)."""
    rows = [list(AP_HEADER)]
    for hz in _bands():
        vals = [0.0] * (len(AP_HEADER) - 1)
        vals[AP_RT30_COL - 1] = (
            float("nan") if hz in STATE["nan_bands"] else STATE["t60"]
        )
        rows.append(["%d Hz" % hz] + vals)
    rows.append(["Global"] + [0.0] * (len(AP_HEADER) - 1))
    rows.append(["Average"] + [0.0] * (len(AP_HEADER) - 1))
    return rows


_INSTALLED = False


def install():
    """Installs the fake modules into sys.modules. Call before importing
    the tools under test.

    Idempotent on purpose: the tools do `import uictrl as ui` at import
    time and hold a reference to that module object. If a second test file
    called install() again it would put a DIFFERENT module in sys.modules
    while the already-imported tools kept pointing at the first one, so
    fixture tweaks would silently not reach the code under test."""
    global _INSTALLED
    if _INSTALLED:
        return sys.modules["uictrl"], sys.modules["libsimpa"]
    _INSTALLED = True
    uictrl = types.ModuleType("uictrl")

    class _ElementType:
        ELEMENT_TYPE_SCENE_SOURCES_SOURCE = "SRC"
        ELEMENT_TYPE_SCENE_RECEPTEURSP = "RECPGRP"
        ELEMENT_TYPE_SCENE_RECEPTEURSP_RECEPTEUR = "RECP"
        ELEMENT_TYPE_REPORT_FOLDER = "REPFOLDER"
        ELEMENT_TYPE_REPORT_GABE_RECP = "REPGABE"
        ELEMENT_TYPE_SCENE_PROJET_ENVIRONNEMENTCONF = "ENVCONF"

    class _IdEvent:
        IDEVENT_RECP_COMPUTE_ACOUSTIC_PARAMETERS = "COMPUTE_AP"
        IDEVENT_RELOAD_FOLDER = "RELOAD"

    uictrl.element_type = _ElementType()
    uictrl.idevent = _IdEvent()
    uictrl._ = lambda s: s

    SCENE_SOURCES = {1: ("SpeakerA", (0.0, 0.0, 2.0))}

    # Each source owns a PROPERTIES child that carries the "enable" bool
    # ("Active source"). The flag is NEVER on the source element itself --
    # Element::GetBoolConfig only sees direct children -- so modelling it
    # on a child is what makes this fixture a real test of the lookup.
    SOURCE_PROPS_ID = {sid: 9000 + sid for sid in SCENE_SOURCES}

    def _scene_receivers():
        # built lazily so tests can change receiver_label after install()
        return {10: (STATE["receiver_label"], (3.0, 0.0, 1.5))}

    class FakeElement:
        def __init__(self, eid):
            self.eid = eid

        def getinfos(self):
            if self.eid == PUNCTUAL_RECEIVERS_GROUP_ID:
                return {"name": "Punctual receivers", "label": "Punctual receivers",
                        "parentid": MODE_ANCESTOR_ID}
            if self.eid == MODE_ANCESTOR_ID:
                return {"name": STATE["mode"], "label": STATE["mode"], "parentid": -1}
            if self.eid == RECEIVER_FOLDER_ID:
                lbl = STATE["receiver_label"]
                return {"name": lbl, "label": lbl, "parentid": PUNCTUAL_RECEIVERS_GROUP_ID}
            if self.eid == SOUND_LEVEL_GABE_ID:
                return {"name": "Sound level", "label": "Sound level",
                        "parentid": RECEIVER_FOLDER_ID,
                        "typeElement":
                            uictrl.element_type.ELEMENT_TYPE_REPORT_GABE_RECP}
            if self.eid in PER_SOURCE_FOLDER_IDS:
                idx = PER_SOURCE_FOLDER_IDS.index(self.eid)
                name = "Source %d" % (idx + 1)
                return {"name": name, "label": name,
                        "parentid": RECEIVER_FOLDER_ID,
                        "typeElement": "REPFOLDER"}
            if self.eid in PER_SOURCE_SOUND_LEVEL_IDS:
                idx = PER_SOURCE_SOUND_LEVEL_IDS.index(self.eid)
                return {"name": "Sound level", "label": "Sound level",
                        "parentid": PER_SOURCE_FOLDER_IDS[idx],
                        "typeElement":
                            uictrl.element_type.ELEMENT_TYPE_REPORT_GABE_RECP}
            if self.eid in SCENE_SOURCES:
                return {"name": SCENE_SOURCES[self.eid][0]}
            if self.eid in SOURCE_PROPS_ID.values():
                return {"name": "Properties", "label": "Properties"}
            if self.eid in _scene_receivers():
                return {"name": _scene_receivers()[self.eid][0]}
            if self.eid == ENVIRONMENT_ID:
                return {"name": "Environment", "label": "Environment"}
            if self.eid == "SCENE_ROOT":
                return {}
            raise KeyError(self.eid)

        def hasproperty(self, name):
            if self.eid == ENVIRONMENT_ID and name == "temperature":
                return True
            # only the PROPERTIES child knows about "enable"
            if name == "enable" and self.eid in SOURCE_PROPS_ID.values():
                return True
            return False

        def getboolconfig(self, name):
            if name == "enable" and self.eid in SOURCE_PROPS_ID.values():
                sid = next(k for k, v in SOURCE_PROPS_ID.items()
                           if v == self.eid)
                return SCENE_SOURCES[sid][0] not in STATE.get(
                    "inactive_sources", ())
            raise KeyError((self.eid, name))

        def getdecimalconfig(self, name):
            if self.eid == ENVIRONMENT_ID and name == "temperature":
                return STATE["temperature_c"]
            raise KeyError((self.eid, name))

        def getpositionconfig(self, key):
            if self.eid in SCENE_SOURCES and key == "pos_source":
                return list(SCENE_SOURCES[self.eid][1])
            if self.eid in _scene_receivers() and key == "pos_recepteur":
                return list(_scene_receivers()[self.eid][1])
            raise KeyError((self.eid, key))

        def childs(self):
            """Direct children only -- this is what distinguishes the
            combined echogram from the per-source ones."""
            if self.eid in SCENE_SOURCES:
                return [(SOURCE_PROPS_ID[self.eid], "F")]
            if self.eid == RECEIVER_FOLDER_ID:
                kids = [(SOUND_LEVEL_GABE_ID, "F")]
                # per-source folders are listed FIRST here on purpose: if a
                # tool ever goes back to a recursive first-match search,
                # this ordering makes it pick the wrong one and fail loudly
                # instead of passing by luck.
                for i in range(_per_source_count()):
                    kids.insert(i, (PER_SOURCE_FOLDER_IDS[i], "F"))
                return kids
            if self.eid in PER_SOURCE_FOLDER_IDS:
                idx = PER_SOURCE_FOLDER_IDS.index(self.eid)
                return [(PER_SOURCE_SOUND_LEVEL_IDS[idx], "F")]
            return []

        def getallelementbytype(self, etype):
            if etype == uictrl.element_type.ELEMENT_TYPE_SCENE_SOURCES_SOURCE:
                return list(SCENE_SOURCES.keys())
            if etype == uictrl.element_type.ELEMENT_TYPE_SCENE_RECEPTEURSP_RECEPTEUR:
                return list(_scene_receivers().keys())
            if (etype == uictrl.element_type.ELEMENT_TYPE_REPORT_GABE_RECP
                    and self.eid == RECEIVER_FOLDER_ID):
                # RECURSIVE, exactly like I-Simpa's own
                # Element::GetAllElementByType (element.cpp:306-316): the
                # per-source reports come back too.
                found = [SOUND_LEVEL_GABE_ID]
                for i in range(_per_source_count()):
                    found.append(PER_SOURCE_SOUND_LEVEL_IDS[i])
                return found
            return []

        def getelementbytype(self, etype):
            if etype == uictrl.element_type.ELEMENT_TYPE_SCENE_RECEPTEURSP:
                return "RECP_GROUP"
            if etype == uictrl.element_type.ELEMENT_TYPE_SCENE_PROJET_ENVIRONNEMENTCONF:
                return ENVIRONMENT_ID if STATE["has_environment"] else -1
            return -1

        def getelementbylibelle(self, name):
            if name == "Acoustic parameters" and self.eid == RECEIVER_FOLDER_ID:
                return ACOUSTIC_PARAMS_ID
            return -1

        # NB: deliberately NO buildfullpath here. In the real API that
        # method exists only on e_file; putting it on element let a real
        # AttributeError slip through to the GUI.

    class FakeEFile(FakeElement):
        """Mirrors uictrl.e_file, which subclasses element and adds
        buildfullpath()."""

        def buildfullpath(self):
            return "C:\\fake\\Punctual receivers\\%s\\" % STATE["receiver_label"]

    def getdataarray(el):
        if el.eid == SOUND_LEVEL_GABE_ID:
            return _make_sound_level_grid()
        if el.eid == ACOUSTIC_PARAMS_ID:
            return _make_acoustic_params_grid()
        raise KeyError(el.eid)

    # NB: attributes are assigned AFTER the class body, because a class
    # body uses LOAD_NAME and therefore cannot see install()'s locals.
    class _Application:
        pass

    _Application.getdataarray = staticmethod(getdataarray)
    _Application.getrootscene = staticmethod(lambda: "SCENE_ROOT")
    _Application.sendevent = staticmethod(lambda el, ev, params=None: None)
    _Application.getuserinput = staticmethod(lambda t, m, f: STATE["user_input"])
    _Application.register_event = staticmethod(lambda fn: "EVT")
    _Application.register_menu_manager = staticmethod(lambda etype, mgr: None)

    uictrl.element = FakeElement
    uictrl.e_file = FakeEFile
    uictrl.application = _Application()
    sys.modules["uictrl"] = uictrl

    libsimpa = types.ModuleType("libsimpa")

    def _require_str(value, where):
        """The real SWIG binding maps std::string to Python 3 str and
        rejects bytes with:
          TypeError: in method 'stringarray_append', argument 2 of type
                     'std::vector< std::string >::value_type const &'
        Passing .encode('cp1252') bytes (as the shipped Python 2-era
        sample scripts do) failed live in the GUI, so the fakes enforce
        the same rule -- otherwise the tests happily accept bytes and the
        bug only shows up in the application."""
        if not isinstance(value, str):
            raise TypeError(
                "in method '%s', argument of type 'std::string' got %s"
                % (where, type(value).__name__)
            )
        return value

    class _StringArray(list):
        def append(self, value):
            list.append(self, _require_str(value, "stringarray_append"))

    class _FloatArray(list):
        def append(self, value):
            if not isinstance(value, float):
                raise TypeError("floatarray_append expects float, got %s"
                                % type(value).__name__)
            list.append(self, value)

    class FakeGabeRw:
        def __init__(self, ncols):
            self.cols = []

        def AppendStrCol(self, arr, name):
            self.cols.append((_require_str(name, "AppendStrCol"), list(arr)))

        def AppendFloatCol(self, arr, name):
            self.cols.append((_require_str(name, "AppendFloatCol"), list(arr)))

        def Save(self, path):
            GABE_WRITES[_require_str(path, "Save")] = self.cols

    libsimpa.Gabe_rw = FakeGabeRw
    libsimpa.stringarray = _StringArray
    libsimpa.floatarray = _FloatArray
    sys.modules["libsimpa"] = libsimpa
    return uictrl, libsimpa


def uictrl_application_getuserinput_swap(new_fn):
    """Swaps the fake getuserinput and returns the previous one, so a test
    can inspect exactly which dialog fields a tool asks for.

    NOTE: assign the plain function, NOT staticmethod(new_fn).
    `uictrl.application` is an INSTANCE, and instance attributes do not go
    through the descriptor protocol -- so a staticmethod object stored here
    is handed back raw when the tool calls it. Python 3.10 made
    staticmethod objects directly callable, but I-Simpa embeds Python
    3.8.1, where that raises "TypeError: 'staticmethod' object is not
    callable". Wrapping it here therefore passes on a modern interpreter
    and fails on the one that actually matters."""
    uictrl = sys.modules["uictrl"]
    previous = uictrl.application.getuserinput
    uictrl.application.getuserinput = new_fn
    return previous


def uictrl_element_class():
    """The fake `element` class, for asserting what it does and does not
    expose (e.g. buildfullpath must live on e_file only)."""
    return sys.modules["uictrl"].element


def written_columns():
    """Returns {row_label: value} for the single receiver column of the
    most recent Gabe write."""
    assert len(GABE_WRITES) == 1, GABE_WRITES
    _, cols = list(GABE_WRITES.items())[0]
    decoded = {}
    for name, values in cols:
        key = name.decode("cp1252") if isinstance(name, bytes) else name
        decoded[key] = [
            v.decode("cp1252") if isinstance(v, bytes) else v for v in values
        ]
    labels = decoded["label"]
    data = decoded[STATE["receiver_label"]]
    return dict(zip(labels, data))
