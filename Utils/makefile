
include ../Compilers.mk

#
# If svnversion utility is found, extract
# SVN revision of working directory. Otherwise,
# insert rev = UNKNOWN.
#
SVN_BIN = svnversion
SVN_LOC = $(shell which $(SVN_BIN))
ifeq (,$(findstring $(SVN_BIN), $(SVN_LOC)))
	SVN_CMD = echo -n "UNKNOWN"
else
	SVN_CMD = echo -n $(strip $(shell $(SVN_BIN) -n .))
endif

BIN_LIST = SVNRevision.py

all: $(BIN_LIST)

SVNRevision.py:
	echo "#!/usr/bin/env python" > SVNRevision.py
	echo "" >> SVNRevision.py
	echo -n 'SVN_VERSION = "' >> SVNRevision.py
	$(SVN_CMD) >> SVNRevision.py
	echo '"'   >> SVNRevision.py

clean :
	-rm -f $(BIN_LIST) SVNRevision.pyc  *~

