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
# TEMPORARY diagnostic tool -- prints the raw getdataarray() structure for
# a receiver's "Sound level" and "Acoustic parameters" reports to the
# Python console, so alcons_tool/sti_tool's parsing assumptions can be
# checked against ground truth instead of guessed at again. Safe to delete
# once alcons_tool/sti_tool are confirmed working.
from __future__ import print_function
import os
import sys
import uictrl as ui

# Make this package folder importable so sibling modules (surface_dump)
# resolve regardless of how I-Simpa's loader set up sys.path.
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)
from libsimpa import *

PUNCTUAL_RECEIVERS_FOLDER_NAME = u"Punctual receivers"


def _print_grid(title, grid, max_rows=6):
    print("----- %s -----" % title)
    if grid is None:
        print("  (None)")
        return
    print("  %d rows total" % len(grid))
    for i, row in enumerate(grid):
        if i >= max_rows:
            print("  ... (%d more rows)" % (len(grid) - max_rows))
            break
        print("  row[%d] (%d items) = %r" % (i, len(row), row))
    if len(grid) > max_rows:
        print("  last row = %r" % (grid[-1],))


def DumpReceiverData(receiver_folder_id):
    receiver = ui.element(receiver_folder_id)
    infos = receiver.getinfos()
    recv_name = infos.get("label", infos.get("name", "?"))
    print("=== Debug dump for receiver '%s' (id=%s) ===" % (recv_name, receiver_folder_id))

    gabe_ids = receiver.getallelementbytype(ui.element_type.ELEMENT_TYPE_REPORT_GABE_RECP)
    print("Children of type ELEMENT_TYPE_REPORT_GABE_RECP: %s" % gabe_ids)
    for gid in gabe_ids:
        child = ui.element(gid)
        print("  id=%s name=%r" % (gid, child.getinfos().get("name")))

    sound_level_id = -1
    for gid in gabe_ids:
        if ui.element(gid).getinfos().get("name") == "Sound level":
            sound_level_id = gid
            break

    if sound_level_id == -1:
        print("No child named exactly 'Sound level' found.")
    else:
        grid = ui.application.getdataarray(ui.element(sound_level_id))
        _print_grid("getdataarray('Sound level')", grid)

        # Trigger RT-30 computation, same call alcons_tool/sti_tool make,
        # then dump 'Acoustic parameters'.
        ui.application.sendevent(
            ui.element(sound_level_id), ui.idevent.IDEVENT_RECP_COMPUTE_ACOUSTIC_PARAMETERS, {"TR": "30"}
        )
        ap_id = receiver.getelementbylibelle("Acoustic parameters")
        print("getelementbylibelle('Acoustic parameters') = %s" % ap_id)
        if ap_id != -1:
            ap_grid = ui.application.getdataarray(ui.element(ap_id))
            _print_grid("getdataarray('Acoustic parameters')", ap_grid)

    print("=== end dump ===")


class manager:
    def __init__(self):
        self.dump_id = ui.application.register_event(self.OnDump)

    def getmenu(self, typeel, idel, menu):
        el = ui.element(idel)
        infos = el.getinfos()
        parentid = infos.get("parentid", -1)
        if parentid is None or parentid < 0:
            return False
        parent = ui.element(parentid)
        if parent.getinfos().get("name") == PUNCTUAL_RECEIVERS_FOLDER_NAME:
            menu.insert(0, ())
            menu.insert(0, ("DEBUG: Dump report data", self.dump_id))
            return True
        return False

    def OnDump(self, idel):
        DumpReceiverData(idel)


ui.application.register_menu_manager(ui.element_type.ELEMENT_TYPE_REPORT_FOLDER, manager())
