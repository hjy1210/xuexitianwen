import numpy as np
from astropy.time import Time
from tools.ecliptic import lonofobj, subtractmod360, nearday, jieqi_index, get_jieqi
def test_lonofsun():
    t = Time('2021-3-20 9:37:28')
    lon = lonofobj('sun', t)
    v = subtractmod360(lon, 0)
    assert  v < 1.0e-5 and v > -1.0e-5
    t = Time('2021-3-20 9:37:29')
    v = subtractmod360(lonofobj('sun', t), 0)
    assert  v < 1.0e-5 and v > -1.0e-5

def test_lonofmoon():
    t = Time('2021-3-20 9:37:28')
    v = subtractmod360(lonofobj('moon', t), 76.6720637)
    assert  v < 1.0e-4 and v > -1.0e-4
    t = Time('2021-3-20 9:37:29')
    v = subtractmod360(lonofobj('moon', t), 76.6722014)
    assert  v < 1.0e-4 and v > -1.0e-4

def test_nearday():
    t = Time('2021-1-25 11:59:59')
    today = nearday(t)
    assert today == Time('2021-1-25')
    t = Time('2021-1-25 12:00:00.1')
    today = nearday(t)
    assert today == Time('2021-1-26')

def test_jieqi_index():
    t = Time('2021-3-20 9:37')
    assert jieqi_index(t) == 23
    t = Time('2021-3-20 9:38')
    assert jieqi_index(t) == 0

def test_get_jieqi():
    assert get_jieqi(Time('2021-3-20 1:1'))['index'] == 0
    assert abs(get_jieqi(Time('2021-3-20 1:1'))['time'].jd - Time('2021-03-20 09:37:28').jd)< 1/86400
    assert get_jieqi(Time('2021-3-19 1:1'))['index'] == -1
    assert get_jieqi(Time('2021-3-19 1:1'))['time'] == None
