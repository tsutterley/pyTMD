import os
import sys
import logging
import subprocess
from setuptools import setup, find_packages

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
log = logging.getLogger()

# get long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

# get install requirements
with open('requirements.txt') as fh:
    install_requires = fh.read().splitlines()

# get version
with open('version.txt') as fh:
    version = fh.read()

# list of all scripts to be included with package
scripts=[os.path.join('scripts',f) for f in os.listdir('scripts') if f.endswith('.py')]

# run cmd from the command line
def check_output(cmd):
    return subprocess.check_output(cmd).decode('utf')

# check if GDAL is installed
gdal_output = [None] * 4
try:
    for i, flag in enumerate(("--cflags", "--libs", "--datadir", "--version")):
        gdal_output[i] = check_output(['gdal-config', flag]).strip()
except:
    log.warning('Failed to get options via gdal-config')
else:
    log.info("GDAL version from via gdal-config: {0}".format(gdal_output[3]))
# if setting GDAL version from via gdal-config
if gdal_output[3]:
    # add version information to gdal in install_requires
    gdal_index = install_requires.index('gdal')
    install_requires[gdal_index] = 'gdal=={0}'.format(gdal_output[3])

setup(
    name='pyTMD',
    version=version,
    description='Tide Model Driver to read OTIS, GOT and FES formatted tidal solutions and make tidal predictions',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/tsutterley/pyTMD',
    author='Tyler Sutterley',
    author_email='tsutterl@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    keywords='Ocean Tides, Load Tides, Pole Tides, Tidal Prediction, OTIS, GOT, FES',
    packages=find_packages(),
    install_requires=install_requires,
    dependency_links=['https://github.com/tsutterley/read-ICESat-2/tarball/main',
        'https://github.com/tsutterley/read-ATM1b-QFIT-binary/tarball/main'],
    scripts=scripts,
    include_package_data=True,
)
