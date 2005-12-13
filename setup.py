from distutils.core import setup, Extension
import os,sys

module1 = Extension('snpindex',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0'),
                                     ('USE_PROJECT_HEADER', '1')],
                    sources = ['snpindex.c', 'cgraph.c'])


if os.access('snpindex.c',os.R_OK):
    print 'Using existing pyrexc-generated C-code...'
else:  # HMM, NO PYREXC COMPILED CODE, HAVE TO RUN PYREXC
    exit_status=os.system('pyrexc snpindex.pyx') # TRY USING PYREX TO COMPILE EXTENSIONS
    if exit_status!=0:  # CAN'T RUN THE PYREX COMPILER TO PRODUCE C
        print '\n\nPyrex compilation failed!  Is pyrex missing or not in your PATH?'
        sys.exit(1) # EXIT WITH ERROR CODE

setup (name = 'snpindex',
       version = '0.1',
       long_description = '''
Common Ancestor Profile functionality as python module.
''',
       ext_modules = [module1])

