from astroquery.jplhorizons import Horizons
from astropy.time import Time
from astropy.coordinates import solar_system_ephemeris, EarthLocation, GeocentricTrueEcliptic, FK5, AltAz, \
    GCRS, CIRS, PrecessedGeocentric, ITRS, TETE
from astropy.coordinates import get_body_barycentric, get_body, get_moon, SkyCoord
import astropy.units as u

import numpy as np

from skyfield.api import load
from skyfield.framelib import ecliptic_frame
from skyfield.api import N, W, E, load, wgs84
ts = load.timescale()
planets = load('de430.bsp')
earth = planets['earth']
moon = planets['moon']
sun = planets['sun']


def test_astroby_skyfield_horizons():
    """
    ### Astropy 與 Skyfield 取得太陽座標對照表

    * 基本上，Astropy 的 get_obj 會得到從地心看太陽的 GCRS 赤道座標，對應於 Skyfield 的 earth.at(t).observe(sun).apparent()。

    * Astropy 用 gcrscoord.transform_to(GeocentricTrueEcliptic(equinox = t)) 將 GCRS 赤道座標轉成真時黃道座標(true ecliptic)，
        對應於 Skyfield 的 gcrscoord.frame_latlon(ecliptic_frame)。
    """
    timestr = "2021-4-20 9:37"
    t_astropy = Time(timestr, scale='utc')
    solar_system_ephemeris.set('de430.bsp')
    sun_astropy_apparent = get_body('sun', t_astropy)
    t_skyfield = ts.utc(2021, 4, 20, 9, 37)
    sun_skyfield_astrometric = earth.at(t_skyfield).observe(sun).radec()
    sun_skyfield_apparent = earth.at(t_skyfield).observe(sun).apparent()
    # print(sun_skyfield_astrometric)
    print("\n比較 GCRS astrometric")
    print("sun_skyfield_astrometric", sun_skyfield_astrometric[0]._degrees,
          sun_skyfield_astrometric[1]._degrees, sun_skyfield_astrometric[2].km)
    # 比較 GCRS
    print("\n比較 GCRS apparent")
    print(sun_astropy_apparent)
    print('sun_astropy_apparent in (deg,deg,km) :', sun_astropy_apparent.ra.deg,
          sun_astropy_apparent.dec.deg, sun_astropy_apparent.distance.km)
    radec_skyfield = sun_skyfield_apparent.radec()
    print('sun_skyfield_apparent in (deg,deg,km) :',
          radec_skyfield[0]._degrees, radec_skyfield[1]._degrees, radec_skyfield[2].km)

    # 比較 GeocentricTrueEcliptic
    print('\n比較 GeocentricTrueEcliptic')
    sun_ecliptic_astropy = sun_astropy_apparent.transform_to(
        GeocentricTrueEcliptic(equinox=t_astropy))
    sun_ecliptic_skyfield = sun_skyfield_apparent.frame_latlon(ecliptic_frame)
    print(sun_ecliptic_astropy)
    print(sun_ecliptic_skyfield[1]._degrees,
          sun_ecliptic_skyfield[0]._degrees, sun_ecliptic_skyfield[2].km)

    # 比較 TETE
    print('\n比較 TreeEclipticTrueEquator(TETE)')
    sun_tete_astropy = sun_astropy_apparent.transform_to(
        TETE(obstime=t_astropy))
    sun_tete_skyfield = sun_skyfield_apparent.radec(epoch='date')
    print(sun_tete_astropy)
    print(sun_tete_skyfield[0]._degrees,
          sun_tete_skyfield[1]._degrees, sun_tete_skyfield[2].km)

    sun_by_earth = Horizons(id='10', id_type='majorbody', location='500', epochs={
                            'start': "2021-4-20 9:37", 'stop': "2021-4-20 9:38", 'step': '1'})
    sun_by_earth_ephem = sun_by_earth.ephemerides(quantities='1,2,4,20,31')
    print("jpl 資料")
    print(sun_by_earth_ephem['RA', 'DEC', 'RA_app',
                             'DEC_app', 'delta', 'ObsEclLon', 'ObsEclLat'])

    # 比較 AltAz
    print("\n比較 AltAz")
    location_astropy = EarthLocation(
        lon=121 * u.deg, lat=25 * u.deg, height=0 * u.m)
    moon_taipei_astropy = get_body('moon', t_astropy, location_astropy)
    # print(moon_taipei_astropy)
    _altaz_astropy = AltAz(obstime=t_astropy, location=location_astropy)
    moon_taipei_altaz_astropy = moon_taipei_astropy.transform_to(
        _altaz_astropy)
    print("moon_taipei_altaz_astropy: ", moon_taipei_altaz_astropy)

    taipei = earth + wgs84.latlon(25 * N, 121 * E, elevation_m=0)
    moon_taipei_skyfield = taipei.at(t_skyfield).observe(moon).apparent()
    alt, az, distance = moon_taipei_skyfield.altaz()
    print("moon_taipei_altaz_skyfield: ",
          az._degrees, alt._degrees,  distance.km)

    moon_by_taipei = Horizons(id='301', id_type='majorbody', location={'lon': 121, 'lat': 25, 'elevation': 0}, epochs={
                              'start': "2021-4-20 9:37", 'stop': "2021-4-20 9:38", 'step': '1'})
    moon_by_taipei_ephem = moon_by_taipei.ephemerides(quantities='1,2,4,20,31')
    print("moon_taipei_jpl (AZ EL)", moon_by_taipei_ephem['AZ', 'EL'])


if __name__ == "__main__":
    test_astroby_skyfield_horizons()
