#!/usr/bin/env python

import os

# Release Tag format: Y.M.Rev
#
# where Y is the year, 9 for 2009, 10 for 2010
#       M is the month 1-12
#       Rev is the mid-month revision beginning with 0 and
#           and incremented by 1
RELEASE_TAG = "9.2.0"
SVN_REV_MODULE = "SVNRevision"


# Empty python class is useful for grouping data like 'C' structs
class Rev:
    pass

try:
    # Attempt import of revision compiled by snvversion
    svn_rev = __import__(SVN_REV_MODULE)
except:
    # Import failed, list version as unknown
    svn_rev = Rev()
    svn_rev.SVN_VERSION = "SVN REV UNKNOWN"

print RELEASE_TAG + " [" + svn_rev.SVN_VERSION + "]"
