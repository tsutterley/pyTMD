import pyTMD
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# adjust style
facecolor = "#fcfcfc"
plt.rcParams["figure.facecolor"] = facecolor
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Lato"]

# angular distnaces from Farrell (1972)
theta = np.array(
    [
        0.0001,
        0.001,
        0.01,
        0.02,
        0.03,
        0.04,
        0.06,
        0.08,
        0.1,
        0.16,
        0.2,
        0.25,
        0.3,
        0.4,
        0.5,
        0.6,
        0.8,
        1.0,
        1.2,
        1.6,
        2.0,
        2.5,
        3.0,
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
        9.0,
        10.0,
        12.0,
        16.0,
        20.0,
        25.0,
        30.0,
        40.0,
        50.0,
        60.0,
        70.0,
        80.0,
        90.0,
        100.0,
        110.0,
        120.0,
        130.0,
        140.0,
        150.0,
        160.0,
        170.0,
        180.0,
    ]
)

# Earth parameters and physical constants
a_axis = pyTMD.earth._wgs84.a_axis
flat = pyTMD.earth._wgs84.flat
GM = pyTMD.earth._wgs84.GM
# universal gravitational constant [N*m^2/kg^2]
G = 6.67430e-11
# average radius of the Earth with same volume as ellipsoid [m]
rad_e = a_axis * np.power(1.0 - flat, 1.0 / 3.0)
# average mass of the Earth
m_e = GM / G
# scale to convert mass to uplift
scale = rad_e / m_e

# Load Love number dataset
lln = "wang-prem"
# calculate scaled Green's functions
GF = scale * pyTMD.earth.greens_function(np.radians(theta), lln=lln)
# normalization in Table A3 of Farrell
norm = 1e12 * rad_e * np.radians(theta)

# create figure
fig, ax = plt.subplots(figsize=(8, 4))
# plot normalized Green's Functions
ax.plot(theta, norm * GF, color="0.4")
ax.set_xlabel("Angular Distance [\u00b0]")
ax.set_ylabel("Green's Function (normalized)")
ax.set_xlim(-1, 181)
ax.xaxis.set_major_locator(MultipleLocator(20))
ax.xaxis.set_minor_locator(MultipleLocator(10))
fig.tight_layout()
plt.show()
