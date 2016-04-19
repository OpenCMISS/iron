# The iron_c.dll/iron.dll files are in the path where iron got installed to;
# we need to add this to the PATH in order to have the import directive find the dependent dlls.
import os
os.environ['PATH'] = r'@NATIVE_LIBRARY_PATH@' + os.pathsep + os.environ['PATH']