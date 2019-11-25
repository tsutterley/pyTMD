from setuptools import setup, find_packages
setup(
	name='Tide Model Driver to read OTIS and GOT formatted tidal solutions and make tidal predictions',
	version='1.0.0.9',
	description='Python',
	url='https://github.com/tsutterley/pyTMD',
	author='Tyler Sutterley',
	author_email='tsutterl@uw.edu',
	license='MIT',
	classifiers=[
		'Development Status :: 3 - Alpha',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Physics',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python :: 2',
		'Programming Language :: Python :: 2.7',
	],
	keywords='Ocean Tides, Load Tides, Pole Tides, Tidal Prediction, OTIS, GOT',
	packages=find_packages(),
	install_requires=['numpy','scipy','pyproj','h5py','netCDF4','matplotlib'],
)
