from distutils.core import setup
from os.path import dirname, join
import imp

with open(join(dirname(__file__), 'README.rst')) as f:
    long_description = f.read()

fn = join(dirname(__file__), 'pyblast.py')
with open(fn) as f:
    pyblast = imp.load_source('pyblast', fn, f)
    version = pyblast.VERSION

setup(
    name='PyBlast',
    version=version,
    description='Run NCBI BLAST with an easy-to-use Pythonic wrapper',
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
