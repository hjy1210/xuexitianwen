
# de430.bsp < de421.bsp < de440.bsp, for de440.bsp difference is VERY large
# 為什麼使用 de440.bsp 時， skyfiled 的計與 jpl 差距如此大？

# The Planetary and Lunar Ephemerides DE430 and DE431
# DE430 and DE431 differ in their integrated time span and lunar dynamical modeling. 
# The dynamical model for DE430 included a damping term between the Moon’s liquid core 
# and solid man-tle that gives the best fit to lunar laser ranging data but that is not 
# suitable for backward integration of more than a few centuries. The ephemeris DE431 is 
# similar to DE430 but was fit without the core/mantle damping term, so the lunar orbit 
# is less accurate than in DE430 for times near the current epoch, but is more suitable 
# for times more than a few centuries in the past. DE431 is a longer integration 
# (covering years –13,200 to +17,191) than DE430 (covering years 1550 to 2650)

# The planetary and lunar ephemeris DE430 succeeds the ephemeris DE421 and its 
# pre-cursor DE405 as a general purpose ephemeris


from skyfield.api import load
from astroquery.jplhorizons import Horizons
from astropy.time import Time
import astropy.units as u
import numpy as np



def compare_obj(planets, code, year, month, day, hour=0, minute=0, second=0):
    '''
    compare position between jpl query and skyfield computation

    parameters:
    ----------
    code : 1,2,3,4,5,6,7,8,9,10,199,299,399,301

    year: int,
    month: int,
    day: int,
    hour: int default 0,
    minute: int default 0,
    second float default 0,
    '''

    # print(code)
    ts = load.timescale()
    t = ts.utc(year, month, day, hour, minute, second)
    objname = planets.names()[code][0]
    # print(objname)
    obj = planets[objname]
    icrf_sf = obj.at(t)
    # print(icrf_sf.position)
    icrf_sf_np = icrf_sf.position.km

    time = Time(t.utc_iso())
    times = time.tt+np.linspace(0, 1, 2)*u.s
    objquery = Horizons(id=str(code), id_type='majorbody', location='@0',
                        epochs={'start': times.value[0], 'stop': times.value[1], 'step': '1'})
    icrf_jpl = objquery.vectors(refplane='earth')
    icrf_jpl['x'].convert_unit_to('km')
    icrf_jpl['y'].convert_unit_to('km')
    icrf_jpl['z'].convert_unit_to('km')
    icrf_jpl_np = np.array(
        [icrf_jpl['x'][0], icrf_jpl['y'][0], icrf_jpl['z'][0]])
    print(objname + ' icrf_sf_np in km', icrf_sf_np)
    print(objname + ' icrf_jpl_np in km', icrf_jpl_np)
    print(objname + ' icrf_sf_np-icrf_jpl_np in km', icrf_sf_np-icrf_jpl_np)
    print()

def testDE(defile):
    print(defile)
    planets = load(defile)
    codes = [code for code in planets.codes if code != 0]
    ts = load.timescale()
    for code in codes:
        compare_obj(planets, code, 2021, 1, 30, 3, 6, 0)


#testDE('de440.bsp')
testDE('de430.bsp')
#testDE('de421.bsp')