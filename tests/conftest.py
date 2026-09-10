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
"""Puts the tool folders on sys.path so the tests can import them.

The tools are installed by copying UserScript/<tool> into I-Simpa's own
UserScript folder, so they are plain directories rather than an installed
package. sti_math.py is kept inside sti_tool/ as the single canonical copy.
"""
import os
import sys

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_USERSCRIPT = os.path.join(_REPO, "UserScript")

for _path in (_USERSCRIPT,
              os.path.join(_USERSCRIPT, "sti_tool"),
              os.path.dirname(os.path.abspath(__file__))):
    if _path not in sys.path:
        sys.path.insert(0, _path)
