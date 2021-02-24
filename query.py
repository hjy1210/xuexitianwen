import astropy.units as u
from astroquery.jplhorizons import Horizons
obj = Horizons(id='301',id_type='majorbody', location=None, epochs = {'start':'2021-3-20 9:36', 'stop':'2021-3-20 9:38', 'step':'1m'})
data = obj.ephemerides(quantities='1,2,4,20,31')
print(type(obj), type(data))
data2 = data['RA','DEC','RA_app', 'DEC_app', 'delta']
print(data2['RA'].quantity, data2['RA'].tolist(),data2['RA'].quantity.value)