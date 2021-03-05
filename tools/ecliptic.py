import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, get_body, GeocentricTrueEcliptic, solar_system_ephemeris
import math
import time


def lonofobj(objname, t):
    # t = Time({'year':year, 'month': month, 'day': day, 'hour': hour, 'minute': minute, 'second': second})
    with solar_system_ephemeris.set('jpl'):
        obj = get_body(objname, t)
    obj_ecliptic = obj.transform_to(GeocentricTrueEcliptic(equinox=t))
    return obj_ecliptic.lon.value


def jieqi_index(t):
    """
    index(0~23) of jiequi time interval
    """
    lon = lonofobj('sun', t)
    return np.floor(lon/15)


def moon_index(t):
    lonMoon = lonofobj('moon', t)
    lonSun = lonofobj('sun', t)
    dif = ((lonMoon - lonSun)+720) % 360
    return np.floor(dif/90)


def get_jieqi_time_in_day(left, right):
    times = left + np.linspace(0, 1, 1441)*(right-left)
    jieqi_indexes = jieqi_index(times)
    for i in range(0, 1440):
        if jieqi_indexes[i+1] != jieqi_indexes[0]:
            times = times[i] + np.linspace(0, 1, 61)*(times[i+1]-times[i])
            jieqi_indexes = jieqi_index(times)
            for j in range(0, 60):
                if jieqi_indexes[j+1] != jieqi_indexes[0]:
                    return {'index': jieqi_indexes[j+1], 'time': times[j]}


def get_suowang_time_in_day(left, right):
    times = left + np.linspace(0, 1, 1441)*(right-left)
    suowang_indexes = moon_index(times)
    for i in range(0, 1440):
        if suowang_indexes[i+1] != suowang_indexes[0]:
            times = times[i] + np.linspace(0, 1, 61)*(times[i+1]-times[i])
            suowang_indexes = moon_index(times)
            for j in range(0, 60):
                if suowang_indexes[j+1] != suowang_indexes[0]:
                    return {'index': suowang_indexes[j+1], 'time': times[j]}


def nearday(t):
    tmp = Time(round(t.jd-0.5)+0.5, format='jd')
    tmp.format = 'iso'
    return tmp


def get_events(start, end, five=True, jieqi=True, suowang=True):
    start = nearday(start)
    end = nearday(end)
    jdlen = end.jd - start.jd
    count = round(jdlen)
    times = Time(start.jd+np.linspace(0, 1, count+1)*jdlen, format='jd')
    jieqis = []
    if jieqi:
        jieqi_indexes = jieqi_index(times)
        for i in range(0, len(jieqi_indexes)-1):
            if jieqi_indexes[i] != jieqi_indexes[i+1]:
                if five:
                    jieqis.append(get_event_time_in_day(
                        times[i], times[i+1], jieqi_index))
                else:
                    jieqis.append(get_jieqi_time_in_day(times[i], times[i+1]))
    suowangs = []
    if suowang:
        suowang_indexes = moon_index(times)
        for i in range(0, len(suowang_indexes)-1):
            if suowang_indexes[i] != suowang_indexes[i+1]:
                if five:
                    suowangs.append(get_event_time_in_day(
                        times[i], times[i+1], moon_index))
                else:
                    suowangs.append(
                        get_suowang_time_in_day(times[i], times[i+1]))

    return jieqis, suowangs


def locate_change_index(indexes):
    left = 0
    right = len(indexes)-1
    for index in range(0, len(indexes)-1):
        if indexes[index+1] != indexes[0]:
            return index


def get_event_time_in_day(left, right, indexmethod):
    times = left + np.linspace(0, 1, 13)*(right-left)
    indexes = indexmethod(times)
    index = locate_change_index(indexes)
    times = times[index] + np.linspace(0, 1, 13)*(times[index+1]-times[index])
    indexes = indexmethod(times)
    index = locate_change_index(indexes)
    times = times[index] + np.linspace(0, 1, 13)*(times[index+1]-times[index])
    indexes = indexmethod(times)
    index = locate_change_index(indexes)
    times = times[index] + np.linspace(0, 1, 11)*(times[index+1]-times[index])
    indexes = indexmethod(times)
    index = locate_change_index(indexes)
    times = times[index] + np.linspace(0, 1, 6)*(times[index+1]-times[index])
    indexes = indexmethod(times)
    index = locate_change_index(indexes)
    return {'index': indexes[index+1], 'time': times[index]}


def subtractmod360(a, b):
    dif = ((a - b) + 720) % 360
    return dif if dif <= 180 else dif - 360


def days_between(start, end):
    start = nearday(start)
    end = nearday(end)
    days = (end-start)/(1.0*u.day)
    count = round(days.value)
    ts = start + np.linspace(0., 1., count + 1)*(end-start)
    return [nearday(t) for t in ts]


def nextday(t):
    today = nearday(t)
    tomorrow = Time(today.jd+1, format='jd')
    tomorrow.format = 'iso'
    return tomorrow


def get_jieqi_time(left, right):
    left_jieqi = jieqi_index(left)
    right_jieqi = jieqi_index(right)
    while (right-left > 1.0*u.second):
        middle = left + 0.5*(right-left)
        middle_jiequi = jieqi_index(middle)
        if middle_jiequi == left_jieqi:
            left = middle
        else:
            right = middle
    return {'index': right_jieqi, 'time': left + 0.5*(right-left)}


def get_moon_time(left, right):
    left_moon = moon_index(left)
    right_moon = moon_index(right)
    while (right-left > 1.0*u.second):
        middle = left + 0.5*(right-left)
        middle_moon = moon_index(middle)
        if middle_moon == left_moon:
            left = middle
        else:
            right = middle
    return {'index': right_moon, 'time': left + 0.5*(right-left)}


def get_jieqi(t):
    """
    get jieqi index(int) and time(Time) in the day of t
    {'index':-1, 'time': None} if not jiequi day
    """
    left = nearday(t)
    right = nextday(left)
    left_jieqi = jieqi_index(left)
    right_jieqi = jieqi_index(right)
    if left_jieqi != right_jieqi:
        return get_jieqi_time(left, right)
    else:
        return {'index': -1, 'time': None}


def get_moon(t):
    """
    get moon index(int) and time(Time) in the day of t
    {'index':-1, 'time': None} if not moon phase day
    """
    left = nearday(t)
    right = nextday(left)
    left_moon = moon_index(left)
    right_moon = moon_index(right)
    if left_moon != right_moon:
        return get_moon_time(left, right)
    else:
        return {'index': -1, 'time': None}


def get_jieqi_between(start, end):
    days = days_between(start, end)
    result = []
    for d in days:
        record = get_jieqi(d)
        if record['index'] >= 0:
            result.append(record)
    return result


def get_moon_between(start, end):
    days = days_between(start, end)
    result = []
    for d in days:
        record = get_moon(d)
        if record['index'] >= 0:
            result.append(record)
    return result


def test1(five=True):
    events = get_events(Time('2020-1-1'), Time('2022-1-1'), five)
    for e in events[0]:
        e['time'].format = 'iso'
        print(e)
    for e in events[1]:
        e['time'].format = 'iso'
        print(e)

def test2():
    result = get_jieqi_between(Time('2020-1-1'), Time('2022-1-1'))  # 驗過
    for r in result:
        print(r)

    result = get_moon_between(Time('2020-1-1'), Time('2022-1-1'))  # 驗過
    for r in result:
        print(r)


if __name__ == "__main__":
    # execute only if run as a script
    import time
    time1 = time.time()
    test1() # 20.941001892089844 seconds
    #test1(False) # 60.31847953796387 seconds
    #test2() # 125.65398097038269 seconds
    time2 = time.time()
    print(time2-time1)
