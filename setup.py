from distutils.core import setup

from beta_rat import __version__

setup(name='beta_rat',
        version=__version__,
        author="Christopher Small",
        author_email="csmall@fhcrc.org",
        scripts=['beta_rat.py'],
        py_modules=['beta_rat'],
        requires=['scipy'])

