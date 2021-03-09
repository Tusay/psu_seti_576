# imports
import numpy as np
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel
from barycorrpy import utc_tdb
from blimpy import Waterfall
from astropy.coordinates import ICRS, SkyCoord
import pdb

def get_drift_rate(file):
    # open the file
    obs = Waterfall(file)

    # get the time of obs + make array 15 min into future
    jdate = obs.header['tstart'] + 2400000.5
    JDUTC = np.linspace(jdate, jdate + (60.0 * 15.0/86400.), num=100)

    # get the pointing of the obs
    c = SkyCoord(obs.header['src_raj'], obs.header['src_dej'])
    s = c.to_string('decimal')
    ra_probe, dec_probe = [float(string) for string in s.split()]

    # other needed params
    obsname = 'GBT'
    epoch = 2451545.0
    rv = 0.0
    zmeas = 0.0
    ephemeris='https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp'

    # get the BC vel
    baryvel = get_BC_vel(JDUTC=JDUTC, ra=ra_probe, dec=dec_probe, obsname='GBT',
                         rv=rv, zmeas=zmeas, epoch=epoch, ephemeris=ephemeris,
                         leap_update=True)

    # take the derivative of velocity to get acceleration (i.e., drift)
    diffT = np.diff(JDUTC) * 86400.0
    diffV = np.diff(baryvel[0])
    drift = diffV/diffT

    # convert m/s^2 to Hz/s and take the max drift calculated
    drift *= (obs.container.f_stop * 1e6 / 3e8)

    # TODO check SIGN!!!
    return np.max(np.abs(drift))
