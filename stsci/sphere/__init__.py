import sys

if sys.version_info[0] >= 3:
    # Python 3 compatibility
    __builtins__['xrange'] = range
