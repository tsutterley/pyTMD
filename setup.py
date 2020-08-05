import os
from setuptools import setup, find_packages

# get long_description from README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

# get install requirements
with open('requirements.txt') as fh:
    install_requires = fh.read().splitlines()

# list of all scripts to be included with package
scripts=[os.path.join('scripts',f) for f in os.listdir('scripts') if f.endswith('.py')]

setup(
    name='pyTMD',
    version='1.0.2.3',
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
        'Programming Language :: Python :: 3.7',
    ],
    keywords='Ocean Tides, Load Tides, Pole Tides, Tidal Prediction, OTIS, GOT, FES',
    packages=find_packages(),
    install_requires=install_requires,
    dependency_links=['https://github.com/tsutterley/read-ICESat-2/tarball/master',
        'https://github.com/tsutterley/read-ATM1b-QFIT-binary/tarball/master'],
    scripts=scripts,
    include_package_data=True,
)
