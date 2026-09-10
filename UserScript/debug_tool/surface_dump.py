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
# Surface-receiver DUMP tool -- groundwork for the STI map.
#
# Adds "DEBUG: Dump surface receiver data" to the right-click menu of
# surface-receiver result folders. It reports the real tree structure and,
# for any surface-receiver result file it can find, the real contents as
# read back through libsimpa's rsurf_io.
#
# WHY A DUMP FIRST: every previous assumption about I-Simpa's data layout
# that was not taken from a live dump turned out wrong (transposed
# echogram axes, 'RT-30' vs 'TR-30', bytes vs str, buildfullpath on the
# wrong class, -inf meaning zero). The STI map needs the per-face echogram
# and a way to write a new map file, so this establishes exactly what is
# there before any of it is built on.
#
# What the STI map will need to confirm from this dump:
#   1. Which tree elements hold surface-receiver results, and their types.
#   2. The on-disk path and extension of the result files.
#   3. record_type, GetTimeStepCount(), GetTimeStep().
#   4. Whether one file holds ONE frequency band (expected) or all of them
#      -- rsurf_data has no frequency dimension in its API, so STI will
#      need one file per octave band, 125 Hz to 8 kHz.
#   5. Per-face record structure: GetFaceRecordCount / GetFaceTimeStep /
#      GetFaceEnergy -- i.e. is a usable per-cell echogram really there.
from __future__ import print_function
import os
import uictrl as ui

try:
    import libsimpa as ls
except ImportError:
    ls = None
    from libsimpa import *  # noqa: F401,F403

SURFACE_RECEIVER_FOLDER_NAMES = (u"Surface receivers", u"Surface receiver")
MAX_LIST = 12


def _type_name_map():
    """Reverse map element_type value -> attribute name, so the dump can
    print readable type names instead of opaque integers."""
    mapping = {}
    et = ui.element_type
    for name in dir(et):
        if name.startswith("ELEMENT_TYPE"):
            try:
                mapping[getattr(et, name)] = name
            except Exception:
                pass
    return mapping


def _property_value(el, infos):
    """Reads a property element's value using the getter matching its type.

    'rstype' in particular decides how the GUI labels a surface map
    (2 = TR -> 's', 9 = STI -> blank, anything unmatched -> 'dB'), so
    seeing its real value is the difference between diagnosing and
    guessing."""
    type_name = ""
    try:
        type_name = str(infos.get("typeElement"))
    except Exception:
        pass
    name = infos.get("name") or u""
    getters = (
        ("getentierconfig", "int"),
        ("getdecimalconfig", "float"),
        ("getboolconfig", "bool"),
        ("getstringconfig", "str"),
    )
    for method, label in getters:
        getter = getattr(el, method, None)
        if getter is None:
            continue
        try:
            parent_id = infos.get("parentid", -1)
            if parent_id is None or parent_id < 0:
                break
            value = getattr(ui.element(parent_id), method)(name)
            return "  value(%s)=%r" % (label, value)
        except Exception:
            continue
    return ""


def _walk(el_id, type_names, depth=0, max_depth=4, seen=None):
    if seen is None:
        seen = set()
    if el_id in seen or depth > max_depth:
        return
    seen.add(el_id)
    el = ui.element(el_id)
    try:
        infos = el.getinfos()
    except Exception as exc:
        print("  " * depth + "id=%s <getinfos failed: %s>" % (el_id, exc))
        return
    tname = type_names.get(infos.get("typeElement"), infos.get("typeElement"))
    print("  " * depth + "id=%-6s name=%-28r type=%-42s%s"
          % (el_id, infos.get("name"), tname, _property_value(el, infos)))
    try:
        children = el.childs()
    except Exception:
        children = []
    for child in children[:MAX_LIST]:
        _walk(child[0], type_names, depth + 1, max_depth, seen)
    if len(children) > MAX_LIST:
        print("  " * (depth + 1) + "... %d more children" % (len(children) - MAX_LIST))


def _candidate_paths(el_id):
    """Every on-disk path we can reach from this element, plus whatever
    files actually live in that folder."""
    paths = []
    try:
        base = ui.e_file(el_id).buildfullpath()
        paths.append(base)
    except Exception as exc:
        print("   e_file(%s).buildfullpath() failed: %s" % (el_id, exc))
        return paths
    folder = base if os.path.isdir(base) else os.path.dirname(base)
    if os.path.isdir(folder):
        try:
            for entry in sorted(os.listdir(folder))[:40]:
                paths.append(os.path.join(folder, entry))
        except Exception as exc:
            print("   listdir failed: %s" % exc)
    return paths


def _try_load_rsurf(path):
    """Attempt to read a surface-receiver result file and dump its real
    structure. Returns True if it loaded."""
    if ls is None:
        print("   libsimpa not importable as a module; cannot test rsurf_io")
        return False
    if not os.path.isfile(path):
        return False
    try:
        data = ls.rsurf_data()
        ok = ls.rsurf_io.Load(path, data)
    except Exception as exc:
        return False
    if not ok:
        return False

    print("   +-- LOADED as rsurf: %s" % path)
    try:
        print("       record type    : %r" % data.getRecordType())
    except Exception as exc:
        print("       record type    : <failed: %s>" % exc)
    try:
        print("       time steps     : %d, step = %g s"
              % (data.GetTimeStepCount(), data.GetTimeStep()))
    except Exception as exc:
        print("       time steps     : <failed: %s>" % exc)
    try:
        print("       nodes          : %d" % data.GetNodesCount())
    except Exception as exc:
        print("       nodes          : <failed: %s>" % exc)

    try:
        rs_count = data.GetRsCount()
    except Exception as exc:
        print("       GetRsCount failed: %s" % exc)
        return True
    print("       surface receivers: %d" % rs_count)

    for rs in range(min(rs_count, 3)):
        try:
            name = data.GetRsName(rs)
            xmlid = data.GetRsXmlId(rs)
            nfaces = data.GetRsFaceCount(rs)
        except Exception as exc:
            print("       rs[%d] <failed: %s>" % (rs, exc))
            continue
        print("       rs[%d] name=%r xmlid=%s faces=%d" % (rs, name, xmlid, nfaces))

        # per-face records: this is the per-cell echogram the STI map needs
        for face in range(min(nfaces, 3)):
            try:
                nrec = data.GetFaceRecordCount(rs, face)
                area = data.ComputeFaceArea(rs, face)
                total = data.GetFaceSumEnergy(rs, face)
            except Exception as exc:
                print("         face[%d] <failed: %s>" % (face, exc))
                continue
            print("         face[%d] records=%d area=%.4f m2 sumEnergy=%g"
                  % (face, nrec, area, total))
            samples = []
            for rec in range(min(nrec, 8)):
                try:
                    step = data.GetFaceTimeStep(rs, face, rec)
                    energy = data.GetFaceEnergy(rs, face, rec)
                    samples.append("(step=%d t=%.4fs E=%g)"
                                   % (step, step * data.GetTimeStep(), energy))
                except Exception as exc:
                    samples.append("<rec %d failed: %s>" % (rec, exc))
            if samples:
                print("            first records: " + " ".join(samples))
    return True


def DumpSurfaceReceivers(folder_id):
    type_names = _type_name_map()
    el = ui.element(folder_id)
    infos = el.getinfos()
    print("=== Surface receiver dump: %r (id=%s) ===" % (infos.get("name"), folder_id))

    print("--- tree structure ---")
    _walk(folder_id, type_names)

    print("--- candidate files, and which load as rsurf ---")
    to_probe = [folder_id]
    try:
        for child in el.childs()[:MAX_LIST]:
            to_probe.append(child[0])
            try:
                for grandchild in ui.element(child[0]).childs()[:MAX_LIST]:
                    to_probe.append(grandchild[0])
            except Exception:
                pass
    except Exception:
        pass

    loaded_any = False
    probed = set()
    for el_id in to_probe:
        for path in _candidate_paths(el_id):
            if path in probed:
                continue
            probed.add(path)
            if os.path.isfile(path):
                marker = "   file: %s (%d bytes)" % (path, os.path.getsize(path))
                if _try_load_rsurf(path):
                    loaded_any = True
                else:
                    print(marker + "   [not an rsurf file]")
    if not loaded_any:
        print("   No file loaded as rsurf. Paths probed:")
        for path in sorted(probed)[:40]:
            print("     %s" % path)
    print("=== end surface dump ===")


class manager:
    def __init__(self):
        self.dump_id = ui.application.register_event(self.OnDump)

    def getmenu(self, typeel, idel, menu):
        el = ui.element(idel)
        infos = el.getinfos()
        name = infos.get("name") or u""
        parent_name = u""
        parentid = infos.get("parentid", -1)
        if parentid is not None and parentid >= 0:
            try:
                parent_name = ui.element(parentid).getinfos().get("name") or u""
            except Exception:
                pass
        if (name in SURFACE_RECEIVER_FOLDER_NAMES
                or parent_name in SURFACE_RECEIVER_FOLDER_NAMES):
            menu.insert(0, ())
            menu.insert(0, ("DEBUG: Dump surface receiver data", self.dump_id))
            return True
        return False

    def OnDump(self, idel):
        DumpSurfaceReceivers(idel)


_mgr = manager()
ui.application.register_menu_manager(ui.element_type.ELEMENT_TYPE_REPORT_FOLDER, _mgr)
try:
    ui.application.register_menu_manager(
        ui.element_type.ELEMENT_TYPE_REPORT_RECEPTEURSSVISUALISATION, _mgr)
except Exception:
    pass
