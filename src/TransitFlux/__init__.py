

from .TransitFlux import *

import sys
if not sys.platform == 'win32':
  from .TransitFlux_ctypes_wrapper import *
