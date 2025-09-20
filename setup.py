from setuptools import Extension, setup
import numpy

module = Extension("_symnmf", sources=['symnmfmodule.c', 'symnmf.c'])
setup(name='_symnmf',
     version='1.0',
     description='Python wrapper for custom C extension',
     ext_modules=[module])