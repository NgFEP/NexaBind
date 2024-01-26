
from __future__ import absolute_import

import os, os.path
import sys
from . import version

if sys.platform == 'win32':
    _path = os.environ['PATH']
    os.environ['PATH'] = '%(lib)s;%(lib)s\plugins;%(path)s' % {
        'lib': version.NexaBind_library_path, 'path': _path}
    try:
        with os.add_dll_directory(version.NexaBind_library_path):
            from . import _NexaBind
    except:
        pass

from NexaBind import *
from NexaBind.vec3 import Vec3



class NexaBindException(Exception):
    """This is the class used for all exceptions thrown by the C++ library."""
    pass
