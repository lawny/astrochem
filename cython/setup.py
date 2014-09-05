# setup.py file
import sys
import os
import shutil

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# clean previous build
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if (name.startswith("astrochem_cython") and not(name.endswith(".pyx") or name.endswith(".pxd"))):
            os.remove(os.path.join(root, name))
    for name in dirs:
        if (name == "build"):
            shutil.rmtree(name)

# build "astrochem_cython.so" python extension to be added to "PYTHONPATH" afterwards...

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension("astrochem",sources=["astrochem_cython.pyx"],
                  libraries=["astrochem"],          # refers to "libexternlib.so"
             )
        ]
)           


