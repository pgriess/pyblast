from distutils.core import setup
from os.path import dirname, join
import imp

# Create the long_description from README.rst, minus the first 6 lines which
# are redundant with the name, version, and short description alreayd displayed
# by PyPi by default.
with open(join(dirname(__file__), 'README.rst')) as f:
    long_description = ''.join(f.readlines()[6:])

fn = join(dirname(__file__), 'pyblast.py')
with open(fn) as f:
    pyblast = imp.load_source('pyblast', fn, f)
    version = pyblast.VERSION

setup(
    name='pyblast',
    version=version,
    description='Run NCBI BLAST with an easy-to-use Pythonic API',
    long_description=long_description,
    author='Peter Griess',
    author_email='pg@std.in',
    py_modules=['pyblast'],
    url='https://github.com/pgriess/pyblast',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries',
    ],
    license='MIT')
