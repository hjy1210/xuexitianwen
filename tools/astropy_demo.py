import astropy.units as u
from astropy.time import Time

def testLeapSecond():
    '''
    2016-12-31 23:59:60 是潤秒，所以 2016-12-31 到 2017-1-1 之間有86401秒
    '''

    print('2016-12-31 23:59:60 是潤秒')
    t1 = Time('2016-12-31')
    t2 = Time('2017-1-1')
    print('t1',t1,'\nt2',t2)
    print('t1.jd:',t1.jd, 't2.jd:', t2.jd, 't2.jd-t1.jd:', t2.jd-t1.jd)
    tdif = t2-t1
    print('((t2-t1).jd-1)*86400: t2-t1 比一天多出的秒數:',(tdif.jd-1)*86400)
    t3 = t1 + 1.*u.day
    print('(t1 + 1.*u.day).ymdhms: ',t3.ymdhms)
    print('# t1.jd, t1.tt.jd, t1.tdb.jd')
    print(t1.jd, t1.tt.jd, t1.tdb.jd)
    print('# t1.jd, (t1.tt.jd-t1.jd)*86400, (t1.tdb.jd-t1.tt.jd)*86400')
    print(t1.jd, (t1.tt.jd-t1.jd)*86400, (t1.tdb.jd-t1.tt.jd)*86400)

if __name__ == "__main__":
    testLeapSecond()