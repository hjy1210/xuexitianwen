
# Astropy vs Horizons
# 利用 astroquery.jplhorizons.Horizons，從 Horizons 網站取得的資料，與 get_body 算得的資料比對。
# RA_app/DEC_app 對應於 get_body 所得的 GCRS 轉成 TETE 座標
# ObsEclLon/ObsEclLat 對應於 get_body 所得的 GCRS 轉成 GeocentricTrueEcliptic 座標

import astropy.units as u
from astropy.time import Time, TimeDelta
from astropy.coordinates import solar_system_ephemeris, EarthLocation, GeocentricTrueEcliptic, FK5, PrecessedGeocentric, TETE, TEME
from astropy.coordinates import get_body_barycentric, get_body, get_moon
from astroquery.jplhorizons import Horizons

def test(timestrs):
    '''
    Astropy vs Horizons
    利用 astroquery.jplhorizons.Horizons，從 Horizons 網站取得的資料，與 get_body 算得的資料比對。
    RA_app/DEC_app 對應於 get_body 所得的 GCRS 轉成 TETE 座標
    ObsEclLon/ObsEclLat 對應於 get_body 所得的 GCRS 轉成 GeocentricTrueEcliptic 座標
    Parameters
    ----------
    timestrs: [iso]
    '''

    # t = Time(["2021-3-20 9:37:27", "2021-3-20 9:37:28", "2021-3-20 9:37:29", "2021-3-20 9:37:30"])
    t = Time(timestrs)
    solar_system_ephemeris.set('de430.bsp')
    # 10 is the code of sun in jpl
    sun = get_body('sun', t)
    sun_geocentrictrueecliptic = sun.transform_to(GeocentricTrueEcliptic(equinox = t))
    sun_tete = sun.transform_to(TETE(obstime = t))
    # 10 is the code of sun in jpl
    sun_by_earth = Horizons(id='10',id_type='majorbody', location='500', epochs = {'start':timestrs[0], 'stop':timestrs[-1], 'step':str(len(timestrs)-1)})
    sun_by_earth_ephem = sun_by_earth.ephemerides(quantities='1,2,4,20,31')
    # column.convert_unit_to is in-place convert
    sun_by_earth_ephem['delta'].convert_unit_to('km')
    print('jpl\n',sun_by_earth_ephem['RA_app','DEC_app','delta','ObsEclLon','ObsEclLat'])
    print('get_body->tete\n', sun_tete)
    print('get_body->geocentrictrueecliptic\n', sun_geocentrictrueecliptic)



if __name__ == "__main__":
    test(["2021-4-20 9:37:27", "2021-4-20 9:37:28", "2021-4-20 9:37:29", "2021-4-20 9:37:30"])