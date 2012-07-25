#    
# Copyright (C) 2005-2011 Hervé Rouault <rouault@lps.ens.fr>
# 
# This file is part of Genherite.
# 
# Genherite is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Genherite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with Genherite; see the file COPYING  If not, see
# <http://www.gnu.org/licenses/>.
#
## Process this file with automake to produce Makefile.in

dist_noinst_DATA = parameters.ggo

cmdline.c: parameters.ggo
	  gengetopt --input=$<

bin_PROGRAMS = genherite

genherite_SOURCES = cmdline.c cmdline.h\
					genherite.cpp genherite.hpp\
				   	evolution.cpp evolution.hpp\
				   	cellule.cpp cellule.hpp\
				   	integr.cpp integr.hpp\
				   	composants.cpp composants.hpp
