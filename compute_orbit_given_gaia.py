import numpy as np
import astropy.units as u
from astropy.table import Table
from galpy.orbit import Orbit
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from galpy.util.conversion import get_physical
from galpy.actionAngle import actionAngleStaeckel, estimateDeltaStaeckel

from galpy.potential.mwpotentials import McMillan17 as pot

#gaiatablefile = '/home/amalie/Nextcloud/projects/roguehennes/gaia/gaiatable_hennes.ascii'
gaiatablefile = '/home/amalie/Nextcloud/projects/maddystar/gaiatable_mygaiaquery.ascii'
a = Table.read(gaiatablefile, format='ascii.commented_header')

useparallax = True

racol = 'RA_GAIA'
deccol = 'DEC_GAIA'
pmracol = 'PMRA_GAIA'
pmdeccol = 'PMDEC_GAIA'
loscol = 'LOS_VELOCITY_GAIA'

if useparallax:
    parallaxcol = 'ZPCORR_PARALLAX_GAIA'
    #parallaxcol = 'PARALLAX_GAIA'
    assert parallaxcol in a.columns
else:
    distcol = 'DISTANCE_BAILERJONES'
    assert distcol in a.columns

for col in [racol, deccol, pmracol, pmdeccol, loscol]:
    assert col in a.columns

ra = a[racol] * u.deg
dec = a[deccol] * u.deg
pmra = a[pmracol] * u.mas / u.year
pmdec = a[pmdeccol] * u.mas / u.year
rv = a[loscol] * u.km / u.second
if useparallax:
    distance = coord.Distance(parallax=a[parallaxcol] * u.mas, allow_negative=True)
else:
    distance = a[distcol] * u.parsec

R0 = get_physical(pot)['ro'] * u.kpc
V0 = get_physical(pot)['vo'] * u.km / u.s

star = SkyCoord(
    ra=ra,
    dec=dec,
    distance=distance,
    pm_ra_cosdec=pmra,
    pm_dec=pmdec,
    radial_velocity=rv,
)

o = Orbit(star, **get_physical(pot))

V = o.V(use_physical=True)
U = o.U(use_physical=True)
W = o.W(use_physical=True)
VPERP = np.sqrt(U**2 + W**2)

print(f"U = {U} km/s, V = {V} km/s, W = {W} km/s, VPERP = {VPERP} km/s")

# Coordinates + velocities in cartesian, Heliocentric frame
star_HC = star.transform_to(coord.Galactic)
star_HC.representation_type = "cartesian"
x_HC = star_HC.u.to(u.kpc)
y_HC = star_HC.v.to(u.kpc)
z_HC = star_HC.w.to(u.kpc)
u_HC = star_HC.U.to(u.km / u.s)
v_HC = star_HC.V.to(u.km / u.s)
w_HC = star_HC.W.to(u.km / u.s)

# Convert to Galactocentric frame
star = star.transform_to(coord.Galactocentric)

# Save coordinates in Galactocentric cartesian coordinates
# The sign of x_GC is changed in order to be in a left-handed frame.
x_GC = -star.x.to(u.kpc)
y_GC = star.y.to(u.kpc)
z_GC = star.z.to(u.kpc)

# Convert to cylindrical Galactocentric coordinates
star.set_representation_cls(
    coord.CylindricalRepresentation, s=coord.CylindricalDifferential
)
R = star.rho.to(u.kpc)
phi = star.phi.to(u.rad)
z = star.z.to(u.kpc)
vR = star.d_rho.to(u.km / u.s)
# Notice the change of sign here, converting to left-handed
# If this is not done, the sign will be wrong in action+angle.
vphi = -(star.d_phi.to(u.rad / u.s) * star.rho.to(u.km) / (1.0 * u.rad))
vz = star.d_z.to(u.km / u.s)

print(x_HC, y_HC, z_HC, u_HC, v_HC, w_HC)
print(R, phi, z)

times = np.linspace(0, 10, 10) * u.Gyr

o.integrate(times, pot)
g_orb = o.getOrbit()

delta = estimateDeltaStaeckel(pot, g_orb[:,0], g_orb[:, 3])
print(f"Focal distance {delta}")
print(o.zmax())
print(o.e())
print(o.E())

aAS = actionAngleStaeckel(
    pot=pot, delta=delta, c=True  # delta=deltas
)

jR, lz, jz, OR, Ophi, Oz, tR, tphi, tz = aAS.actionsFreqsAngles(
    R, vR, vphi, z, vz, phi, c=True
)

# A proper determination of guiding centre radii would require the
# integration of the stellar orbits within an assummed potential
# However, under the assumption of a flat rotation curve for the
# Galaxy, we can approximate R_guide as
rguide = (lz * R0 * V0) / V0
# Here we assumme the circular velocity to be that of the Sun.

# galpy ouputs the actions and frequencies in its internal coordinates
# where distances in units of R0 and V0 so if we want the quantities
# in physical units, we need to scale them back
action_unit = R0.value * V0.value
freq_unit = (1 / R0.value) * V0.value
jR *= action_unit
lz *= action_unit
jz *= action_unit
OR *= freq_unit
Ophi *= freq_unit
Oz *= freq_unit
tR = tR * u.rad
tphi = tphi * u.rad
tz = tz * u.rad

# Restrict the azimuthal angle to be within 0-2pi.
# tphi[tphi > np.pi] -= 2 * np.pi

# Compute orbital parameters such as eccentricity and zmax from Staeckel
# approximation
e, zmax, rperi, rap = aAS.EccZmaxRperiRap(R, vR, vphi, z, vz, phi, use_physical=True)

print(jR, lz, jz, rguide, e, zmax)
