try:
    from setuptools import setup, Extension
except ImportError:
    print 'Setuptools not imported, falling back to distutils'
    from distutils.core import setup, Extension

import os,sys

if 'setuptools' in sys.modules:
    cmdclass = {}
else:
    from Pyrex.Distutils import build_ext
    cmdclass = {'build_ext': build_ext}

module1 = Extension('snpindex.core',
                    define_macros = [('MAJOR_VERSION', '0'),
                                     ('MINOR_VERSION', '1'),
                                     ('USE_PROJECT_HEADER', '1')],
                    sources = [os.path.join('snpindex', 'core.pyx'),
                               os.path.join('snpindex', 'cgraph.c')])



setup (name = 'snpindex',
       version = '0.1',
       long_description = '''Fast SNP indexing and pair analysis.''',
       author = "Christopher Lee",
       author_email='leec@chem.ucla.edu',
       packages = ['snpindex'],
       ext_modules = [module1],
       cmdclass = cmdclass,
       )

