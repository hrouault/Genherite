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

ACLOCAL_AMFLAGS = -I m4

BUILT_SOURCES = $(top_srcdir)/.version

$(top_srcdir)/.version:
	echo $(VERSION) > $@-t && mv $@-t $@
dist-hook:


SUBDIRS = src doc

EXTRA_DIST = autogen.sh

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	@if test -d "$(srcdir)/.git"; \
		then \
		echo Creating ChangeLog && \
		( cd "$(top_srcdir)" && \
		echo '# Generated by Makefile. Do not edit.'; echo; \
		$(top_srcdir)/missing --run git log --stat ) > ChangeLog.tmp \
		&& mv -f ChangeLog.tmp $(distdir)/ChangeLog \
		|| ( rm -f ChangeLog.tmp ; \
		echo Failed to generate ChangeLog >&2 ); \
		else \
		echo A git clone is required to generate a ChangeLog >&2; \
		fi

