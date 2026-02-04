
try:
    import cupy
    np = cupy
except ImportError:
    import numpy
    np = numpy
