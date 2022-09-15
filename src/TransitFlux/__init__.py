
from .TransitFlux import *

if not sys.platform == 'win32':
  from .TransitFlux_ctypes_wrapper import *
