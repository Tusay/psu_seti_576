# imports
import numpy as np
from math import isclose

import skyfield
from skyfield.api import load, Star, T0
from skyfield.data import hipparcos
from skyfield.positionlib import ICRF, position_of_radec

import pytz
from datetime import datetime, timedelta, timezone
from pytz import timezone
from tzlocal import get_localzone

# function to get probe coords opposite hipparcos star
def get_probe_coords(hip_id, distance):
    # set timezone + load planet info
    utcTZ = timezone("UTC")
    local_tz = get_localzone()
    ts = load.timescale()
    planets = load('de421.bsp')

    # get time of observation and locations of observation
    t = ts.from_datetime(datetime(2020, 11, 9, 00, 00,tzinfo=utcTZ))
    sun = planets['sun']
    earth = planets['earth']

    # load in the hipparcos data
    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)

    # # match on star coords
    # ra_inds = []
    # for ind, cat_ra in enumerate(df['ra_degrees'].values):
    #     if isclose(cat_ra, ra_star, rel_tol=1e-2):
    #         ra_inds.append(ind)

    # dec_inds = []
    # for ind, cat_dec in enumerate(df['dec_degrees'].values):
    #     if isclose(cat_dec, dec_star, rel_tol=1e-2):
    #         dec_inds.append(ind)

    # # find where ra & dec match within tolerance
    # sets = set(ra_inds).intersection(set(dec_inds))
    # inds = [x for x in iter(sets)]
    # print(inds)

    # just take the first one for now
    the_star = Star.from_dataframe(df.loc[hip_id])

    # "observe" the star from Earth
    starPosObj = earth.at(t).observe(the_star)
    ra, dec, dist = starPosObj.radec()

    # get coords of opposite point on the sky
    ra = ((ra._degrees + 180.0) % 360) * (24.0/360.0)
    dec = -1.0 * dec._degrees

    # get coords of the probe and return ra dec in degree
    probe = position_of_radec(ra, dec, distance, t=t)
    sunPosObj = earth.at(t).observe(sun)
    inverseSunCoord = -1.0 * sunPosObj.position.au
    probeCoord = probe.position.au
    probeDir = ICRF(probeCoord + inverseSunCoord).radec()
    return probeDir[0]._degrees, probeDir[1]._degrees

# alpha_cen_hip_id = 71683

# ra, dec = get_probe_coords(alpha_cen_hip_id, 50.0)
