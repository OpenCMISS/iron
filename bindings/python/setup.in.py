#!/usr/bin/env python

import os
from sys import platform
from setuptools import setup
from setuptools.dist import Distribution

requires = ['numpy']
package_data = {'opencmiss.iron': [@SETUP_PY_PACKAGE_FILES_STR@]}


class BinaryDistribution(Distribution):
    def is_pure(self):
        return False

    def has_ext_modules(self):
        return True


setup(
    name='opencmiss.iron',
    version='@Iron_VERSION@@IRON_DEVELOPER_VERSION@',
    description=('Python bindings for the OpenCMISS computational '
            'modelling library Iron.'),
    long_description=('Python bindings to OpenCMISS-Iron. '
            'OpenCMISS (Open Continuum Mechanics, Imaging, Signal processing '
            'and System identification) is a mathematical modelling '
            'environment that enables the application of finite element '
            'analysis techniques to a variety of complex '
            'bioengineering problems. OpenCMISS-Iron is the computational backend component '
            'of OpenCMISS.'),
    author='Hugh Sorby',
    license='Mozilla Tri-license',
    author_email='hsorby@auckland.ac.nz',
    url='http://www.opencmiss.org/',
    requires=requires,
    packages=['opencmiss', 'opencmiss.iron'],
    package_data=package_data,
    distclass=BinaryDistribution,
    include_package_data=True,
    zip_safe=False,
)
