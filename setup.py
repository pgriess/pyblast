from distutils.core import setup
from os.path import dirname, join

with open(join(dirname(__file__), 'README.rst')) as f:
    long_description = f.read()

setup(
    name='PyBlast',
    version='0.1',
    description='Run NCBI BLAST with an easy-to-use Pythonic wrapper',
    long_description=long_description,
    author='Peter Griess',
    author_email='pg@std.in',
    py_modules=['pyblast'])
