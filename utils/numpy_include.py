import numpy
import sys

try:
    print(sys.version)
    print(numpy.get_include())
except AttributeError:
    print(numpy.get_numpy_include())
