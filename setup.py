"""
setup.py: Used for building python wrappers for test NexaBind library.
"""
import ast
import re
import os
import sys
import platform
import numpy
#from distutils.core import setup
from setuptools import setup
from Cython.Build import cythonize

MAJOR_VERSION_NUM='1'
MINOR_VERSION_NUM='0'
BUILD_INFO='0'
IS_RELEASED = False




def reportError(message):
    sys.stdout.write("ERROR: ")
    sys.stdout.write(message)
    sys.stdout.write("\nExiting\n")
    sys.exit(1)


def buildKeywordDictionary(major_version_num=MAJOR_VERSION_NUM,
                           minor_version_num=MINOR_VERSION_NUM,
                           build_info=BUILD_INFO):
    from setuptools import Extension
    setupKeywords = {}
    setupKeywords["name"]              = "NexaBind"
    setupKeywords["version"]           = "%s.%s.%s" % (major_version_num,
                                                       minor_version_num,
                                                       build_info)
    setupKeywords["author"]            = "IQB Group"
    setupKeywords["license"]           = \
    "Python Software Foundation License (BSD-like)"
    #setupKeywords["url"]               = "https://NexaBind.org"
    #setupKeywords["download_url"]      = "https://NexaBind.org"
    setupKeywords["packages"]          = [
                                          "NexaBind",
                                          "NexaBind.unit",
                                          "NexaBind",
                                          "NexaBind.app",
                                          "NexaBind.app.internal"]
    setupKeywords["data_files"]        = []
    setupKeywords["package_data"]      = {"NexaBind" : [],
                                          "NexaBind.app" : ['data/*.xml', 'data/*.pdb', 'data/amber14/*.xml', 'data/charmm36/*.xml', 'data/implicit/*.xml'],
                                          "NexaBind.app.internal" : []}
    setupKeywords["platforms"]         = ["Linux", "Mac OS X", "Windows"]
    setupKeywords["description"]       = \
    "Python wrapper for NexaBind (a C++ MD package)"
    setupKeywords["long_description"]  = \
    """NexaBind is a toolkit for molecular simulation. 
    """

    return setupKeywords


def main():
    if sys.version_info < (2, 7):
        reportError("NexaBind requires Python 2.7 or better.")
    if platform.system() == 'Darwin':
        macVersion = [int(x) for x in platform.mac_ver()[0].split('.')]
        if tuple(macVersion) < (10, 5):
            reportError("NexaBind requires Mac OS X Leopard (10.5) or better.")

    setupKeywords=buildKeywordDictionary()
    setup(**setupKeywords)

if __name__ == '__main__':
    main()
