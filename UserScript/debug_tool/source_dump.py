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
# Scene SOURCE dump -- finds out exactly where the "Active source" flag
# lives, instead of inferring it.
#
# WHY: alcons_tool must count only ACTIVE sources. A disabled source is
# never written into the calculation XML (e_scene_sources_source.h:104-111)
# so SPPS never traces it; giving it a direct-sound arrival steals energy
# from E_rev into E_dir and makes %ALcons look far too good.
#
# The C++ says the flag is a BOOL named "enable" on the source's PROPERTIES
# child (e_scene_sources_source_properties.h:49), read via
# GetBoolConfig, which only sees DIRECT children (element.cpp:1300-1314).
# This dump CHECKS that against the live tree rather than trusting it --
# every previous structural assumption in this project that was not dumped
# turned out wrong.
from __future__ import print_function
import uictrl as ui

MAX_CHILD = 30


def _try(label, fn):
    try:
        return "%s=%r" % (label, fn())
    except Exception as exc:
        return "%s=<%s>" % (label, exc)


def _probe(el_id, indent):
    """Everything we can learn about one element's 'enable' flag."""
    pad = "  " * indent
    try:
        el = ui.element(el_id)
    except Exception as exc:
        print(pad + "id=%s <element() failed: %s>" % (el_id, exc))
        return
    try:
        infos = el.getinfos()
    except Exception as exc:
        print(pad + "id=%s <getinfos failed: %s>" % (el_id, exc))
        return
    print(pad + "id=%-6s name=%-22r type=%-40s %s  %s"
          % (el_id, infos.get("name"), infos.get("typeElement"),
             _try("hasproperty('enable')",
                  lambda: el.hasproperty("enable")),
             _try("getboolconfig('enable')",
                  lambda: el.getboolconfig("enable"))))


def DumpSceneSources():
    print("=== SCENE SOURCE dump: where does 'Active source' live? ===")
    try:
        rootscene = ui.element(ui.application.getrootscene())
    except Exception as exc:
        print("   could not reach the scene root: %s" % exc)
        return
    try:
        src_ids = rootscene.getallelementbytype(
            ui.element_type.ELEMENT_TYPE_SCENE_SOURCES_SOURCE)
    except Exception as exc:
        print("   could not list sources: %s" % exc)
        return

    print("   %d source element(s) found in the scene" % len(src_ids))
    for sid in src_ids:
        print("   ---------------------------------------------")
        _probe(sid, 1)
        try:
            children = ui.element(sid).childs()
        except Exception as exc:
            print("      childs() failed: %s" % exc)
            continue
        print("      %d direct child element(s):" % len(children))
        for child in children[:MAX_CHILD]:
            _probe(child[0], 3)
            # one more level: the BOOL itself may sit under Properties
            try:
                grandkids = ui.element(child[0]).childs()
            except Exception:
                grandkids = []
            for gk in grandkids[:MAX_CHILD]:
                _probe(gk[0], 4)
    print("=== end scene source dump ===")


class manager:
    def __init__(self):
        self.dump_id = ui.application.register_event(self.OnDump)

    def getmenu(self, typeel, idel, menu):
        menu.insert(0, ())
        menu.insert(0, ("DEBUG: Dump scene sources", self.dump_id))
        return True

    def OnDump(self, idel):
        DumpSceneSources()


_mgr = manager()
for _type_name in ("ELEMENT_TYPE_SCENE_SOURCES",
                   "ELEMENT_TYPE_SCENE_SOURCES_SOURCE"):
    try:
        ui.application.register_menu_manager(
            getattr(ui.element_type, _type_name), _mgr)
    except Exception:
        pass
