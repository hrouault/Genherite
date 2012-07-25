#    
# Copyright (C) 2006-2012 Hervé Rouault <rouault@lps.ens.fr>
# Copyright (C) 2012 Marc Santolini <santolin@lps.ens.fr>
# 
# This file is part of Genherite.
# 
# Imogene is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Imogene is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Imogene; see the file COPYING  If not,
# see <http://www.gnu.org/licenses/>.

EXTRA_DIST = genherite.pyx genherite.c

genherite.c: genherite.pyx
	$(PYTHON) $(srcdir)/setup.py build_ext --force --inplace

pyexec_LTLIBRARIES = genherite.la
nodist_genherite_la_SOURCES = genherite.c
genherite_la_LDFLAGS = -avoid-version -module
