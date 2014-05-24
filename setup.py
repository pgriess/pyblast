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
    license='MIT')
