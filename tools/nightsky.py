from astropy.io import ascii
import astropy.units as u
from astropy.time import Time
import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.coordinates import Angle, GeocentricTrueEcliptic, GCRS
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, solar_system_ephemeris, get_body #, TETE
from astropy.table import QTable, Table, Column, vstack
from enum import Enum
from matplotlib.font_manager import FontProperties
import datetime

def getPlanetTable(timestr, location):
    objnames = ['sun','moon', 'mercury','venus', 'mars', 'jupiter', 'saturn','uranus', 'neptune']
    Vmag = [-3,-3,-2.48,-4.92,-2.94, -2.94, -0.55, 5.38, 7.67]
    size = [17,17,5,5,5,4,4,2,1]
    t = Time(timestr, scale='utc')
    solar_system_ephemeris.set('jpl')
    table = Table ()
    alts=[]
    azs=[]
    distances = []
    for i in range(len(objnames)):
        position = get_body(objnames[i], t, location=location)
        altaz = position.transform_to(AltAz(obstime=t, location=location))
        alts.append(altaz.alt.degree)
        azs.append(altaz.az.degree)

    table['name']= objnames
    table['alt']=alts
    table['az']=azs
    table['Vmag']=Vmag
    table['size'] = size
    return table

def getEclipticTable(timestr):
    lon=[i*15 for i in range(25)]
    dec=[0 for i in range(25)]
    ecliptic= SkyCoord(lon=np.array(lon)*u.deg, lat=np.array(dec)*u.deg, equinox=timestr, frame='geocentrictrueecliptic')
    gcrs = ecliptic.transform_to(GCRS(obstime=timestr))
    table = Table()
    table["RA"]=gcrs.ra.degree
    table["DE"]=gcrs.dec.degree
    return table
class Projection(Enum):
    Perspective = 1
    Lambert = 2
    Stereo = 3

def vmag2size(vmag):
    size = np.exp((4-vmag)/2*np.log(2.512))
    return size

def localtime2utctime(localtimestr, timezone):
    time=Time(localtimestr)-timezone*u.hour
    time.format='iso'
    return time.value
# theta 方位角，phi 仰角
def spher2cart(theta,phi):
      return np.array([np.cos(phi)*np.sin(theta),np.cos(phi)*np.cos(theta),np.sin(phi)])

def spher2cartDphi(theta,phi):
  return np.array([-np.sin(phi)*np.sin(theta),-np.sin(phi)*np.cos(theta),np.cos(phi)])

def get_project_matrix(azindegree,altindegree):
    altrad = altindegree*np.pi/180
    azrad = azindegree*np.pi/180
    e2 = spher2cart(azrad,altrad)
    e3 = spher2cartDphi(azrad,altrad)
    e1 = np.cross(e2,e3)
    M = np.stack([e1,e2,e3])
    #print(e1,e2,e3)
    #print(M)
    return M

def set_size(table):
    table['size']=vmag2size(table['Vmag'])

bscTable = ascii.read("catalogs\\bsc5\\bsc5.dat", format='fixed_width_no_header',
      names=('name', 'RAh', 'RAm', 'RAs','DE-', 'DEd', 'DEm', 'DEs', 'Vmag', 'pmRA', 'pmDE'),
      col_starts=(14,75,77,79,83,84,86,88,102,148,154),
      col_ends=  (24,76,78,82,83,85,87,89,106,153,159)             
      )
bscTable['RA']=(bscTable['RAh']+bscTable['RAm']/60+bscTable['RAs']/3600)*15
bscTable['DE']=(bscTable['DEd']+bscTable['DEm']/60+bscTable['DEs']/3600)
for i in range(len(bscTable)):
    if bscTable[i]['DE-']=='-':
        bscTable[i]['DE']= -bscTable[i]['DE']
set_size(bscTable)

orionTable = ascii.read("catalogs\\onions.dat")
orionTable['RA'] = Angle(np.array(orionTable['RA'])).degree
orionTable['DE'] = Angle(np.array(orionTable['DE'])).degree

def set_altaz(table, obstime, location):
    #if 'distance' in table.colnames:
    #    stars = SkyCoord(ra=table['RA']*u.deg, dec = table['DE']*u.deg, distance=table['distance']*u.km, frame='tete', obstime=obstime)
    #else:
    #    stars = SkyCoord(ra=table['RA']*u.deg, dec = table['DE']*u.deg, frame='gcrs')
    stars = SkyCoord(ra=table['RA']*u.deg, dec = table['DE']*u.deg, frame='gcrs')
    altaz = AltAz(obstime=obstime, location=location)
    stars_by_taipei = stars.transform_to(altaz)

    table['alt']=stars_by_taipei.alt
    table['az']=stars_by_taipei.az

def set_xyz_proj(table, az0, alt0, projection:Projection, fov):
    cart = spher2cart(np.array(table['az'])*np.pi/180,np.array(table['alt'])*np.pi/180)
    M = get_project_matrix(az0,alt0)
    coords = M.dot(cart)
    table['x']=coords.T[:,0]
    table['y']=coords.T[:,1]
    table['z']=coords.T[:,2]

    hfov=fov/2*np.pi/180
    #view_table = brightest_table[brightest_table['y']> np.cos(hfov)]
    #view_table = table
    ##### view_table = view_table[view_table['alt']>0]
    if projection == Projection.Stereo:
        # stereographic projection
        table['x'] = 2*table['x'] /(1+table['y'])
        table['z'] = 2*table['z'] /(1+table['y'])
        #scale = np.sin(hfov)/(1+np.cos(hfov))/np.sqrt(2)
        scale = 2*np.tan(hfov/2)
    elif projection == Projection.Lambert:
        # Lambert euqal area projection
        r = np.sqrt(table['x']*table['x']+table['z']*table['z'])
        s = np.sqrt(2-2*table['y'])
        table['x'] = table['x'] * s/r
        table['z'] = table['z'] * s/r
        #scale = np.sqrt(1-np.cos(hfov))
        scale = np.sqrt(2-2*np.cos(hfov))
    else:
        # perspective projection
        table['x']= table['x']/table['y']
        table['z']= table['z']/table['y']
        #scale = np.tan(hfov)/np.sqrt(2)
        scale = np.tan(hfov)
    return scale

def proj(table, projection:Projection, Vmag_max, obstime, location, az0, alt0, fov):
    def on_pick(event):
        line = event.artist
        #xdata, ydata = line.get_data()
        ind = event.ind
        print('on pick line:', ind, table['name'][ind[0]])

    #mask = table['Vmag'] < Vmag_max
    #brightest_table = table[mask]
    #brightest_table = table
    set_altaz(table, obstime, location)

    scale = set_xyz_proj(table, az0, alt0, projection, fov)



    fig, ax = plt.subplots(figsize=(16, 9)) 
    ax.scatter(table['x'],table['z'], s=(table['size']*(72./fig.dpi))**2, c='w',marker='.', picker=5)
    ax.set_xlim(-16/9*scale,16/9*scale)
    ax.set_ylim(-1*scale,1*scale)
    ax.set_aspect('equal')
    # using chinese characters
    #prop = FontProperties()
    #prop.set_file('c:\\windows\\fonts\\mingliu.ttc')
    #plt.title('{}, ({}, {}), ({}, {}, {}), {}楊宏章'.format(projection, location.lon, location.lat, az0,alt0,fov, obstime), fontproperties=prop)
    ax.set_title('{}, ({}, {}), ({}, {}, {}), {}'.format(projection, location.lon, location.lat, az0,alt0,fov, obstime))
    ax.set_facecolor((0, 0, 0))
    ax.set_xticks([])
    ax.set_yticks([])
    #cid = fig.canvas.mpl_connect('button_press_event', onclick)
    cid = fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    
def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))


if __name__ == "__main__":
    obstime = '2021-2-16 13:32'
    location = EarthLocation(lon=121*u.deg, lat= 25*u.deg)
    Vmag_max = 6
    az0 = 0     # 0/180
    alt0 = 25    # 25/-25
    fov = 120
    #proj(table, Projection.Perspective, Vmag_max, obstime, location, az0, alt0, fov)
    #proj(table, Projection.Lambert, Vmag_max, obstime, location, az0, alt0, fov)
    #proj(bscTable, Projection.Stereo, Vmag_max, obstime, location, az0, alt0, fov)
    #proj(orionTable, Projection.Perspective, Vmag_max, obstime, location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 13:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 14:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 15:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 16:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 17:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 18:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 19:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 20:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 21:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 22:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-16 23:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 0:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 1:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 2:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 3:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 4:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 5:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 6:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 7:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 8:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 9:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 10:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 11:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 12:32', location, az0, alt0, fov)
    #proj(table, Projection.Stereo, Vmag_max, '2021-2-17 13:32', location, az0, alt0, fov)
    #table = getPlanetTable(obstime)
    #print(table)
    #print(vstack([orionTable,table2]))
    #timestr ='{}'.format(datetime.datetime.now())
    #s = localtime2utctime(timestr,8)
    #print(type(s),s)
    obstime = '2021-02-23 13:57:0'
    location = EarthLocation(lon=121.5528*u.deg, lat=24.99136*u.deg)
    table= getPlanetTable(obstime,location)
    print(table)
    t = Time(obstime)
    print(t)
    #self.lon = '121.5528'
    #self.lat = '24.99136'
    sun = get_body('sun',t, location)
    sun_altaz= sun.transform_to(AltAz(obstime=t, location=location))
    print(sun_altaz)

