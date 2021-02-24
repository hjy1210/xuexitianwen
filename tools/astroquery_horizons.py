import numpy as np
import astropy.units as u
from astropy.constants import c
from astropy.time import Time
from astropy.coordinates import SkyCoord, TETE, Angle
from astroquery.jplhorizons import Horizons


def test_Ephemerides_Vectors_Astrometric():
    strstartdate = '2021-1-1 12:00:00'
    strenddate = '2021-1-2 12:00:00'
    start_time = Time(strstartdate)
    start_time_tdb = start_time.tdb
    end_time = Time(strenddate)
    end_time_tdb = end_time.tdb
    # print(start_time.value, start_time_tdb.value)
    mars_by_earth_obs = Horizons(id='499', id_type='majorbody', location='500',
                                 epochs={'start': start_time.value, 'stop': end_time.value, 'step': '1'})
    mars_by_earth_vec = Horizons(id='499', id_type='majorbody', location='500',
                                 epochs={'start': start_time_tdb.value, 'stop': end_time_tdb.value, 'step': '1'})
    mars_by_earth_obs_data = mars_by_earth_obs.ephemerides(
        quantities='1,2,4,20,31')
    mars_by_earth_vec_astrometric = mars_by_earth_vec.vectors(
        refplane='earth', aberrations='astrometric', delta_T=True)
    mars_by_earth_vec_apparent = mars_by_earth_vec.vectors(
        refplane='earth', aberrations='apparent', delta_T=True)
    mars_by_earth_vec_geometric = mars_by_earth_vec.vectors(
        refplane='earth', delta_T=True)
    print("比較 Horizons 中，用 Observer 模式查得的 RA/DEC與用Vectors模式(aberrations='astrometric')查得的是否相同？是相同！")
    print(mars_by_earth_obs_data['datetime_str','RA', 'DEC', 'RA_app', 'DEC_app', 'delta'])
    # print(mars_by_earth_vec_astrometric['x','y','z','lighttime'])
    astrometric_position = SkyCoord(x=mars_by_earth_vec_astrometric['x'].tolist()*u.au,
                                    y=mars_by_earth_vec_astrometric['y'].tolist()*u.au,
        z=mars_by_earth_vec_astrometric['z'].tolist()*u.au, representation_type='cartesian', frame='gcrs')
    astrometric_position.representation_type = 'spherical'
    print(astrometric_position)

    # print("compare TETE")
    # apparent_position = SkyCoord(x = mars_by_earth_vec_apparent['x'].tolist()*u.au,
    #    y = mars_by_earth_vec_apparent['y'].tolist()*u.au,
    #    z = mars_by_earth_vec_apparent['z'].tolist()*u.au, representation_type='cartesian', frame='gcrs')
    # apparent_position_tete = apparent_position.transform_to(TETE(obstime = start_time.value))
    # apparent_position.representation_type='spherical'    
    # print(apparent_position)
    # apparent_position_tete.representation_type='spherical'
    # print(apparent_position_tete)
    # print(mars_by_earth_vec_geometric['datetime_jd','x','y','z','lighttime'])

def test_Ephemerides_Vectors_Astrometric2():
    print("比較 Horizons 中，用 Observer 模式查得的 RA/DEC與用Vectors模式(aberrations='geometric')查得的加上光線旅行時間是否相同？相同！")
    strstartdate = '2021-1-1 12:00:00'
    strenddate = '2021-1-2 12:00:00'
    start_time = Time(strstartdate)
    start_time_tdb = start_time.tdb
    end_time = Time(strenddate)
    end_time_tdb = end_time.tdb
    # print(start_time.value, start_time_tdb.value)
    mars_by_earth_obs = Horizons(id='499', id_type='majorbody', location='500',
                                 epochs={'start': start_time.value, 'stop': end_time.value, 'step': '1'})

    mars_by_earth_obs_data = mars_by_earth_obs.ephemerides(quantities='1,2,4,20,31')

    print(mars_by_earth_obs_data['datetime_str', 'RA','DEC','RA_app','DEC_app','delta'])
    earth_by_bary_vec = Horizons(id='399', id_type='majorbody', location='@0',
                                 epochs={'start': start_time_tdb.value, 'stop': end_time_tdb.value, 'step': '1'})
    earth_position = earth_by_bary_vec.vectors(refplane='earth')
    earth_current_position=np.array([earth_position['x'][0],earth_position['y'][0],earth_position['z'][0]])
    
    jdlag = ((mars_by_earth_obs_data[0]['delta']*u.au).to(u.m)/c).value/86400
    jdnow = start_time_tdb.jd
    jd0 = jdnow - jdlag

    start_time_tdb = Time(jd0-1/1440, format='jd', scale='tdb')
    start_time_tdb.format = 'iso'
    end_time_tdb = Time(jd0+1/1440, format='jd', scale='tdb')
    end_time_tdb.format = 'iso'


    mars_by_bary_vec = Horizons(id='499', id_type='majorbody', location='@0',
                                 epochs={'start': start_time_tdb.value, 'stop': end_time_tdb.value, 'step': '120'})
    mars_position = mars_by_bary_vec.vectors(refplane='earth')

    found = False
    for i in range(len(mars_position)):
        v = np.array([mars_position['x'][i],mars_position['y'][i],mars_position['z'][i]])
        v = v-earth_current_position
        d = np.linalg.norm(v)
        lighttime = (((d*u.au).to(u.m))/c).value/86400
        # print(i,'jdnow-mars_position["datetime_jd"][i]:', jdnow-mars_position['datetime_jd'][i],' lighttime:',lighttime)
        if jdnow-mars_position['datetime_jd'][i] <= lighttime:
            astrometric = SkyCoord(x=v[0]*u.au, y=v[1]*u.au,
                                       z=v[2]*u.au, representation_type='cartesian', frame='gcrs')
            astrometric.representation_type = 'spherical'
            print(astrometric)
            found = True
            break
    if not found:
        print('search failed')


def deg2hms(d):
    d = d/15
    hour = np.floor(d)
    m = (d-hour)*60
    min = np.floor(m)
    s = (m-min)*60
    return str(hour)+'h'+str(min)+'m'+str(s)+'s'
def deg2dms(d):
    deg = np.floor(d)
    m = (d-deg)*60
    min = np.floor(m)
    s = (m-min)*60
    return str(deg)+'d'+str(min)+'\''+str(s)+'"'


if __name__ == "__main__":
    test_Ephemerides_Vectors_Astrometric()
    test_Ephemerides_Vectors_Astrometric2()
    #print(deg2hms(25.00173), deg2hms(25.0032375))
    #print(deg2dms(11.32678), deg2dms(11.32744477))
    #print(Angle(25.00173, unit=u.deg).hms,Angle(25.0032375, unit=u.deg).hms)
    #print(Angle(11.32678, unit=u.deg).dms, Angle(11.32744477, unit=u.deg).dms)
    #print(Angle(11.32678, unit=u.deg).value)
    #print(np.outer(np.ones(4),[5,6]))
