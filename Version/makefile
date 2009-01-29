
include ../Compilers.mk

BIN_LIST = SVNRevision.py

all: $(BIN_LIST)

SVNRevision.py:
	echo "#!/usr/bin/env python" > SVNRevision.py
	echo "" >> SVNRevision.py
	echo -n 'SVN_VERSION = "' >> SVNRevision.py
	svnversion -n .. >> SVNRevision.py
	echo '"'   >> SVNRevision.py

clean :
	-rm -f $(BIN_LIST)
