from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='pyTMD',
    version='1.0.1.7',
    description='Tide Model Driver to read OTIS and GOT formatted tidal solutions and make tidal predictions',
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
    keywords='Ocean Tides, Load Tides, Pole Tides, Tidal Prediction, OTIS, GOT',
    packages=find_packages(),
    install_requires=['numpy','scipy','pyproj','h5py','netCDF4','matplotlib','cartopy','ipyleaflet'],
    dependency_links=['https://github.com/tsutterley/read-ICESat-2/tarball/master',
        'https://github.com/tsutterley/read-ATM1b-QFIT-binary/tarball/master'],
)
