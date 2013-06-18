from setuptools import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

from os import path

from betarat.version import __version__


ext_modules = [Extension(
    name="rf_hyp2f1",
    sources=[path.join("betarat/rf_hyp2f1", x) for x in ('rf_hyp2f1.pyx', 'rf_hyp2f1_core.c')],
    include_dirs=['.'],
    language="c",
    libraries=['gmp']
    )]

setup(name='betarat',
        version=__version__,
        url='http://github.com/fhcrc/betarat',
        author="Christopher Small",
        author_email="csmall@fhcrc.org",
        entry_points={'console_scripts': ['betarat = betarat.scripts.cli:main']},
        packages=['betarat', 'betarat/rf_hyp2f1', 'betarat/scripts'],
        cmdclass={"build_ext": build_ext},
        ext_modules=ext_modules,
        requires=['scipy', 'mpmath'])

