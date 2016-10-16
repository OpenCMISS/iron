#!/usr/bin/env python

import os
from sys import platform
from setuptools import setup

requires = []#['numpy']
package_data = {'opencmiss.iron': ['$<TARGET_FILE_NAME:@IRON_PYTHON_MODULE@>']}

#try:
    #if platform == 'darwin':
    #    os.symlink('@IRON_TARGET_FILE@', 'opencmiss/iron/libiron.dylib')
    #    package_data['opencmiss.iron'].append('libiron.dylib')

setup(
    name='OpenCMISS-Iron',
    version='@Iron_VERSION@',
    description=('Python bindings for the OpenCMISS computational '
            'modelling library Iron.'),
    long_description=('Python bindings to OpenCMISS-Iron. '
            'OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing '
            'and System identification) is a mathematical modelling '
            'environment that enables the application of finite element '
            'analysis techniques to a variety of complex '
            'bioengineering problems. OpenCMISS-Iron is the computational backend component '
            'of OpenCMISS.'),
    author='Adam Reeve',
    license='Mozilla Tri-license',
    author_email='hsorby@aucklanduni.ac.nz',
    url='http://www.opencmiss.org/',
    install_requires=requires,
    packages=['opencmiss', 'opencmiss.iron'],
    package_data=package_data
)
#finally:
    
#    if platform == 'darwin':
#        os.unlink('opencmiss/iron/libiron.dylib')
