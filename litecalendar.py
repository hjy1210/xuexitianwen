from astropy.time import Time,TimeDelta
from astropy.coordinates import get_sun,FK5,PrecessedGeocentric,get_body,Angle,solar_system_ephemeris, GeocentricTrueEcliptic
from astropy import coordinates as coord
from astropy import constants as const
from astropy.coordinates import SkyCoord
import astropy.units as u
import math
import numpy as np
jieqinames=["春分","清明","穀雨","立夏","小滿","芒種","夏至","小暑","大暑","立秋",
              "處暑","白露","秋分","寒露","霜降","立冬","小雪","大雪","冬至","小寒","大寒","立春","雨水","驚蟄"]
moonshapenames=['朔月','上弦','望月','下弦']
Cnumber=['一','二','三','四','五','六','七','八','九','十','十一','十二']
codedict= {'O': '全在光亮區','H':'部分在半影區', 'F': '部分在本影區', 'G': '全在本影區','P':'部分在半影區',
           'D':'全在偽本影','C':'部分在偽本影'}

def CCTerm(start,end,type):
    eighthours=TimeDelta(8*3600,format='sec')
    def termtime(daytime,daytime2,code,unit='hour'):
        if unit=='hour':
            times=daytime+np.linspace(0.,1.,25)*(daytime2-daytime)
        elif unit=='minute':
            times=daytime+np.linspace(0.,1.,61)*(daytime2-daytime)
        else:
            times=daytime+np.linspace(0.,1.,61)*(daytime2-daytime)
        gmttimes=times-eighthours
        if type=='jieqi':
            diflongitudes=get_body('sun', gmttimes,  ephemeris='jpl').\
                transform_to(GeocentricTrueEcliptic(equinox=gmttimes)).lon.deg
        else:
            sunlongitudes=get_body('sun', gmttimes,  ephemeris='jpl').\
                transform_to(GeocentricTrueEcliptic(equinox=gmttimes)).lon.deg
            moonlongitudes=get_body('moon', gmttimes,  ephemeris='jpl').\
                transform_to(GeocentricTrueEcliptic(equinox=gmttimes)).lon.deg
            diflongitudes=(moonlongitudes-sunlongitudes)%360
        difcode=diflongitudes // sec
        difcode=[int(difcode[i]) for i in range(len(difcode))]
        for i in range(len(difcode)-1):
          if difcode[i+1]!=code:
            if unit=='hour':
              return termtime(times[i],times[i+1],code,"minute")
            else:
              if unit=='minute':
                return termtime(times[i],times[i+1],code,"second")
              else:
                return ((code+1) % (360//sec),(times[i]+0.5*(times[i+1]-times[i])).value)
        return (-1,"")

    start=Time(Time(start,out_subfmt='date').iso)
    end=Time(Time(end,out_subfmt='date').iso)
    daycount=int(round((Time(end)-Time(start)).value))
    # print(start,end,daycount)
    # print(daycount)
    times=Time(start)+np.linspace(0.,1.,daycount+1)*(end-start+TimeDelta(3,format='sec')) # 考量潤秒，故加上3秒
    times=Time(Time(times,out_subfmt='date').iso)
    gmttimes=times-TimeDelta(8*3600,format='sec')
    if type=='jieqi':
        diflongitudes= get_body('sun', gmttimes,  ephemeris='jpl').\
            transform_to(GeocentricTrueEcliptic(equinox=gmttimes)).lon.deg
        sec=15
    else:
        sunlongitudes=get_body('sun', gmttimes,  ephemeris='jpl').\
            transform_to(GeocentricTrueEcliptic(equinox=gmttimes)).lon.deg
        moonlongitudes=get_body('moon', gmttimes,  ephemeris='jpl').\
            transform_to(GeocentricTrueEcliptic(equinox=gmttimes)).lon.deg
        diflongitudes=(moonlongitudes-sunlongitudes)%360
        sec=90
    res=[]
    difcode=diflongitudes // sec
    difcode=[int(difcode[i]) for i in range(len(difcode))]
    for i in range(len(difcode)-1):
      if difcode[i]!=difcode[i+1]:
        res.append(termtime(times[i],times[i+1],difcode[i],'hour'))
    return res

def _CCalendar(year):
    oneday=TimeDelta(1, format='jd')
    limits=CCTerm(str(year-1)+'-12-1',str(year+1)+'-1-1','jieqi')
    limits=[data for data in limits if data[0]==18] # 冬至
    darkmoon=CCTerm(Time(limits[0][1])+oneday,Time(limits[1][1])+oneday,'moon')
    
    darkmoon=[m[1] for m in darkmoon if m[0]==0]  # 朔月
    
    res=[]
    
    if len(darkmoon)<13:
        for i in range(12):
            j=(i+12)%12
            if j==0:
                j=12
            res.append((j,darkmoon[i]))
    else:
        add=False
        mnumber=12
        for i in range(12):
            js=CCTerm(Time(darkmoon[i]),Time(darkmoon[i+1]),'jieqi')
            #print(js)
            js=[d for d in js if d[0]%2==0] # 中氣
            #print(js)
            if add==False and len(js)==0 :
                res.append((mnumber-1+0.1,darkmoon[i])) # 有小數點代表閏月
                add=True
            else:
                res.append((mnumber,darkmoon[i]))
                mnumber=mnumber+1
                if mnumber==13:
                    mnumber=1
        res.append((mnumber,darkmoon[12]))
    return(res)
def CCalendar(year):
    res=_CCalendar(year)+_CCalendar(year+1)
    return([m for m in res if m[1]>=str(year)+'-01-01' and m[1]<str(year+1)+'-01-01'])

def PrintCCalendar(year):
    res=CCalendar(year)
    for m in res:
        name=Cnumber[int(m[0])-1]+'月初一'
        if int(m[0])!=m[0]:
            name='潤'+name
        print(name+','+m[1])
    
# http://eclipse.gsfc.nasa.gov/eclipse.html 可查到歷年的日蝕月蝕的資料
oneminute=TimeDelta(60, format='sec')
eighthours=TimeDelta(8*3600,format='sec')
onesecond=TimeDelta(1, format='sec')
def separation(a,b):
    """a,b兩向量的夾角"""
    val=a.dot(b)/(np.linalg.norm(a)*np.linalg.norm(b))
    #if val>1.0:
    #    val=1.0
    #if val< -1.0:
    #    val= -1.0
    return(math.acos(val))
def ShadowCone(s,rs,o,ro):
    o2s=s-o
    #dso=np.array([np.linalg.norm(x) for x in o2s])
    dso=np.linalg.norm(o2s,axis=1)
    sodir=np.array([-o2s[i]/dso[i] for i in range(len(o2s))])
    #sodir=np.diag(1.0/dso).dot(-o2s) ##### inefficient
    du=ro/(rs-ro)*dso
    ucenter=np.array([o[i]+du[i]*sodir[i] for i in range(len(o))])
    #ucenter=o+np.diag(du).dot(sodir) ##### inefficient
    #uarc=np.array([math.asin(rs/(dso[i]+du[i])) for i in range(len(s))])
    #if asin:
    uarc=np.arcsin(ro/du)
    #else:
    #    uarc=ro/du
    dp=ro/(rs+ro)*dso
    #parc=np.array([math.asin(rs/(dso[i]-dp[i])) for i in range(len(s))])
    #if asin:
    parc=np.arcsin(ro/dp)
    #else:
    #    parc=ro/dp
    pcenter=np.array([o[i]-dp[i]*sodir[i] for i in range(len(o))])
    #pcenter=o-np.diag(dp).dot(sodir) ##### inefficient
    return((ucenter,pcenter,uarc,parc,sodir))
def Region(S,O,T,rt,Ucenter,Pcenter,uarc,parc,SOdir):
    res=[]
    for i in range(len(T)):
        if (T[i]-O[i]).dot(SOdir[i])<=0:
            res.append('O')    # O 代表光亮區
            continue
        vp=T[i]-Pcenter[i]
        ptarc=math.asin(rt/np.linalg.norm(vp))   # 從 P 看標的球圓盤T的視角之半, delta
        psep=separation(T[i]-Pcenter[i],SOdir[i]) # \angle TPO
        if psep>parc[i]+ptarc: # TPO>beta+delta
            res.append('O')   # psep>parc[i]+ptarc 代表標的球完全落在光亮區
            continue
        vu=T[i]-Ucenter[i]
        dist=np.linalg.norm(vu)
        if vu.dot(SOdir[i])<=0:  # T 在 U 的右側
            if dist<=rt:         # 本影頂點在球T內部
                res.append('F')  # 球T與本影區有重疊
                continue
            utarc=math.asin(rt/dist) # gamma
            usep=separation(vu,-SOdir[i]) #usep=angle TUO,alpha=uarc[i]
            if usep<= uarc[i]-utarc: # TUO<=alpha-gamma
                res.append('G')  # 球T完全進入本影區
            elif usep<= uarc[i]+utarc: # alpha-gamma<TUO<=alpha+gamma
                res.append('F')  # 球T部分落在本影區
            else:               # alpha+gamma < TUO
                res.append('P')  # 球T與半影區有重疊
            continue
        else:                    # T 在 U 的左側
            if dist<=rt:         # 本影頂點在球T內部
                res.append('F')  # 球T與本影區有重疊
                continue
            utarc=math.asin(rt/dist) #gamma
            usep=separation(vu,SOdir[i]) #pi-TUO
            if usep<= uarc[i]-utarc:  # pi-TUO < alpha-gamma , TUO>pi-alpha+gamma
                res.append('D')  # 球T全部落在偽本影
            elif usep<= uarc[i]+utarc: # alpha-gamma<pi-TUO < alpha+gamma, pi-alpha-gamma< TUO <pi-alpha+gamma
                res.append('C')  # 球T局部落在偽本影
            else:                # pi-TUO>alpha+gamma, TUO<pi-alpha-gamma
                res.append('P')  # 球T與半影區有重疊
            continue
    return(res)  
def EarthRegion(times):
    gmttimes=times-eighthours
    sun=get_body('sun', gmttimes,  ephemeris='jpl')
    moon=get_body('moon', gmttimes,  ephemeris='jpl')
    earth=get_body('earth', gmttimes,  ephemeris='jpl')
    rs=const.R_sun.to('km').value
    re=const.R_earth.to('km').value
    ro=1737.1
    sunpos=np.transpose(np.array([sun.cartesian.x.value,
                                  sun.cartesian.y.value,sun.cartesian.z.value]))
    moonpos=np.transpose(np.array([moon.cartesian.x.value,
                                  moon.cartesian.y.value,moon.cartesian.z.value]))
    earthpos=np.zeros((len(sunpos),3))
    Ucenter,Pcenter,uarc,parc,SOdir=ShadowCone(sunpos,rs,moonpos,ro)
    return(Region(sunpos,moonpos,earthpos,re,Ucenter,Pcenter,uarc,parc,SOdir))
def SolarEclipseOne0(start,end):
    #print((Time(end)-Time(start)).value)
    minutecount=int(np.ceil((Time(end)-Time(start)).value*1440))
    times=Time(start)+np.linspace(0.,1.,minutecount+1)*(Time(end)-Time(start))
    regions=EarthRegion(times)
    #print(times[0],regions[0])
    for i in range(len(times)-1):
        if regions[i]!=regions[i+1]:
            newsymbol=regions[i+1]
            oldsymbol=regions[i]
            sectimes=times[i]+onesecond*np.linspace(0,60,61)
            secregions=EarthRegion(sectimes)
            for j in range(61):
                j1=60-j
                if secregions[j1]!=newsymbol:
                    print(sectimes[j1],newsymbol)
                    break
def SolarEclipseOne(start,end):
    #print((Time(end)-Time(start)).value)
    minutecount=int(np.ceil((Time(end)-Time(start)).value*1440))
    times=Time(start)+np.linspace(0.,1.,minutecount+1)*(Time(end)-Time(start))
    regions=EarthRegion(times)
    #print(times[0],regions[0])
    for i in range(len(times)-1):
        if regions[i]!=regions[i+1]:
            newsymbol=regions[i+1]
            oldsymbol=regions[i]
            sectimes=times[i]+(times[i+1]-times[i])*np.linspace(0.,1.,61)
            secregions=EarthRegion(sectimes)
            for j in range(61):
                j1=60-j
                if secregions[j1]!=newsymbol:
                    print(sectimes[j1],'地球',codedict[oldsymbol],'->',codedict[newsymbol])
                    break
def SolarEclipse(start,end):
    darkmoon=CCTerm(start,end,'moon')
    darkmoon=[codedate[1] for codedate in darkmoon if codedate[0]==0]
    #print(darkmoon)
    for i in range(len(darkmoon)):
        halfday=TimeDelta(0.5,format="jd")
        dark=Time(darkmoon[i])
        start=dark-halfday
        end=dark+halfday
        #print([str(start),str(end)])
        SolarEclipseOne(str(start),str(end))
def LunarRegion(times,air):
    gmttimes=times-eighthours
    sun=get_body('sun', gmttimes,  ephemeris='jpl')
    moon=get_body('moon', gmttimes,  ephemeris='jpl')
    earth=get_body('earth', gmttimes,  ephemeris='jpl')
    rs=const.R_sun.to('km').value
    re=const.R_earth.to('km').value
    rm=1737.1
    sunpos=np.transpose(np.array([sun.cartesian.x.value,
                                  sun.cartesian.y.value,sun.cartesian.z.value]))
    moonpos=np.transpose(np.array([moon.cartesian.x.value,
                                  moon.cartesian.y.value,moon.cartesian.z.value]))
    earthpos=np.zeros((len(sunpos),3))
    # reference http://eclipse.gsfc.nasa.gov/LEcat5/shadow.html
    # NASA 用 Danjon 的方法，將地球放大成 1.01倍，中央氣象局用 Chauvenet 的方法，類似於將地球放大1.02倍
    # 此處採用Danjon 的方法
    Ucenter,Pcenter,uarc,parc,SOdir=ShadowCone(sunpos,rs,earthpos,(1+air)*re)
    return(Region(sunpos,earthpos,moonpos,rm,Ucenter,Pcenter,uarc,parc,SOdir))
def LunarEclipseOne(start,end,air=0.02):
    minutecount=int(np.ceil((Time(end)-Time(start)).value*1440))
    times=Time(start)+np.linspace(0.,1.,minutecount+1)*(Time(end)-Time(start))
    regions=LunarRegion(times,air)
    # print(times[0],regions[0])
    for i in range(len(times)-1):
        if regions[i]!=regions[i+1]:
            oldsymbol=regions[i]
            newsymbol=regions[i+1]
            sectimes=times[i]+onesecond*np.linspace(0,60,61)
            secregions=LunarRegion(sectimes,air)
            for j in range(61):
                j1=60-j
                if secregions[j1]!=newsymbol:
                    print(sectimes[j1],'月亮',codedict[oldsymbol],'->',codedict[newsymbol])
                    break
def LunarEclipse(start,end,air=0.02):
    darkmoon=CCTerm(start,end,'moon')
    darkmoon=[codedate[1] for codedate in darkmoon if codedate[0]==2]
    #print(darkmoon)
    for i in range(len(darkmoon)):
        halfday=TimeDelta(0.5,format="jd")
        dark=Time(darkmoon[i])
        start=dark-halfday
        end=dark+halfday
        #print([str(start),str(end)])
        LunarEclipseOne(str(start),str(end))
def Eclipse(start,end,air=0.02):
    darkmoon=CCTerm(start,end,'moon')
    # darkmoon=[codedate[1] for codedate in darkmoon if codedate[0]==0 ]
    #print(darkmoon)
    for i in range(len(darkmoon)):
        if darkmoon[i][0] % 2 ==1 : continue
        halfday=TimeDelta(0.5,format="jd")
        dark=Time(darkmoon[i][1])
        start=dark-halfday
        end=dark+halfday
        #print([str(start),str(end)])
        if darkmoon[i][0]==0:
            SolarEclipseOne(str(start),str(end))
        else:
            LunarEclipseOne(str(start),str(end),air)
def obj_sight0(date,objname,diskanglestr,lonstr,latstr,tzinhours):
    """
    date: string, date of observing. 
    objname : string, name of sun, moon or planet in solar system.
    diskanglestr : string, apparent angle of objname.
    lonstr : string, longitute of observing location on earth.
    latstr : string, latitute of observing location on earth.
    ltzinhours : float, time zone in hours of observing location.
    return : print out rise time, set time, transit time and cooresponding alt/az.
    """
    lon=Angle(lonstr)
    lat=Angle(latstr)
    halfdisk=Angle(diskanglestr)*0.5
    date=Time(Time(date,out_subfmt='date').iso)
    times=date+np.linspace(0.,1440.,1441)*TimeDelta(60,format='sec')  # 每分鐘計算位置
    deltat=tzinhours*TimeDelta(3600,format='sec')
    a=const.R_earth.to('m').value       # 地球半徑
    with solar_system_ephemeris.set('jpl'):
        n=len(times)
        gmttimes=times-deltat   # 世界時間 UTC
        # 標的物在真赤道座標系統座標
        bodies=get_body(objname,gmttimes).transform_to(PrecessedGeocentric(equinox=gmttimes))  
        sides=gmttimes.sidereal_time('apparent', 'greenwich')+lon    # local sidereal time
        x=bodies.cartesian.x.to('m').value        # x,y,z為標的物在及時赤道座標系統的直角座標   
        y=bodies.cartesian.y.to('m').value
        z=bodies.cartesian.z.to('m').value
        px=np.array(a*np.cos(lat)*np.cos(sides))  # px,py,pz 為觀測地在真赤道座標系統的直角坐標
        py=np.array(a*np.cos(lat)*np.sin(sides))
        pz=np.array(a*np.sin(lat)*np.ones(n))
    positions=np.vstack((px,py,pz)).T  # n x 3 array
    objs=np.vstack((x,y,z)).T          # n x 3 array
    # e_z 的每一行為天頂方向單位向量，e_n的每一行為地表往北單位向量，e_e的每一行為地表往東單位向量
    e_z=np.vstack((np.cos(lat)*np.cos(sides),np.cos(lat)*np.sin(sides),np.sin(lat)*np.ones(n)))  # 3 x n
    e_n=np.vstack((-np.sin(lat)*np.cos(sides),-np.sin(lat)*np.sin(sides),np.cos(lat)*np.ones(n))) # 3 x n
    e_e=np.vstack((-np.sin(sides),np.cos(sides),np.zeros(n)))                 # 3 x n
    coords=np.zeros((n,3))   # coords 的每一列代表標的物在地表座標系統的直角坐標
    az=np.zeros(n)           # az 標的物的方位角，北方為0，往東為正，正東方90度，正南方180度，正西方270度。
    alt=np.zeros(n)          # alt 標的物的仰角
    for i in range(n):
        coords[i]=np.array(np.vstack((e_e[:,i],e_n[:,i],e_z[:,i]))).dot(objs[i]-positions[i])
        r=np.sqrt(coords[i][0]**2+coords[i][1]**2)
        rho=np.sqrt(coords[i][0]**2+coords[i][1]**2+coords[i][2]**2)
        alt[i]=np.arcsin(coords[i][2]/rho)*180/np.pi
        az[i]=np.arctan2(coords[i][0],coords[i][1])*180/np.pi
        if az[i]<0:
            az[i]=az[i]+360
    thresh=-(halfdisk+Angle("0d34m0s")).deg   # 大氣折射，讓光線轉彎了34角分
    altdiff=alt-thresh        # altdif 變號的時候就是標的物升起或沉沒的區間，用內插進行估計
    for i in range(len(times)-1):
        if altdiff[i]*altdiff[i+1] <=0 :
            t=times[i]+(times[i+1]-times[i])*(0-altdiff[i])/(altdiff[i+1]-altdiff[i])
            if altdiff[i]<0:
                print("Rise",t,"Az",az[i])
            else:
                print("Set",t,"Az",az[i])
        if not (az[i] // 180 == az[i+1] // 180) and alt[i]>0:
            if np.abs(az[i]-180)<30:
                ratio=(180-az[i])/(az[i+1]-az[i])
                t=times[i]+(times[i+1]-times[i])*ratio
                print("Transit",t,"Alt",alt[i]+(alt[i+1]-alt[i])*ratio,"Az","South")
            else:
                left=az[i]
                right=az[i+1]
                if left>300:
                    left=left-360
                if right>300:
                    right=right-360
                ratio=(0-left)/(az[i+1]-az[i])
                t=times[i]+(times[i+1]-times[i])*ratio
                print("Transit",t,"Alt",alt[i]+(alt[i+1]-alt[i])*ratio,"Az","North")
def obj_sight(date,objname,lonstr,latstr,tzinhours,approx=True):
    """
    date: string, date of observing. 
    objname : string, name of sun, moon or planet in solar system.
    lonstr : string, longitute of observing location on earth.
    latstr : string, latitute of observing location on earth.
    ltzinhours : float, time zone in hours of observing location.
    approx : bool, True for delta_ut1_utc=0
    return : print out rise time, set time, transit time and cooresponding alt/az.
    """
    lon=Angle(lonstr)
    lat=Angle(latstr)
    if objname=='sun':
      r_obj=const.R_sun.to('m').value
    else:
      if objname=="moon":
        r_obj=1737000
    #halfdisk=Angle(diskanglestr)*0.5
    date=Time(Time(date,out_subfmt='date').iso)
    times=date+np.linspace(0.,1440.,1441)*TimeDelta(60,format='sec')  # 每分鐘計算位置
    deltat=tzinhours*TimeDelta(3600,format='sec')
    a=const.R_earth.to('m').value       # 地球半徑
    with solar_system_ephemeris.set('jpl'):
        n=len(times)
        gmttimes=times-deltat   # 世界時間 UTC
        # 標的物在真赤道座標系統座標
        bodies=get_body(objname,gmttimes).transform_to(PrecessedGeocentric(equinox=gmttimes)) 
        if approx:
          gmttimes.delta_ut1_utc = 0. 
        sides=gmttimes.sidereal_time('apparent', 'greenwich')+lon    # local sidereal time
        x=bodies.cartesian.x.to('m').value        # x,y,z為標的物在及時赤道座標系統的直角座標   
        y=bodies.cartesian.y.to('m').value
        z=bodies.cartesian.z.to('m').value
        px=np.array(a*np.cos(lat)*np.cos(sides))  # px,py,pz 為觀測地在真赤道座標系統的直角坐標
        py=np.array(a*np.cos(lat)*np.sin(sides))
        pz=np.array(a*np.sin(lat)*np.ones(n))
    positions=np.vstack((px,py,pz)).T  # n x 3 array
    objs=np.vstack((x,y,z)).T          # n x 3 array
    # e_z 的每一行為天頂方向單位向量，e_n的每一行為地表往北單位向量，e_e的每一行為地表往東單位向量
    e_z=np.vstack((np.cos(lat)*np.cos(sides),np.cos(lat)*np.sin(sides),np.sin(lat)*np.ones(n)))  # 3 x n
    e_n=np.vstack((-np.sin(lat)*np.cos(sides),-np.sin(lat)*np.sin(sides),np.cos(lat)*np.ones(n))) # 3 x n
    e_e=np.vstack((-np.sin(sides),np.cos(sides),np.zeros(n)))                 # 3 x n
    coords=np.zeros((n,3))   # coords 的每一列代表標的物在地表座標系統的直角坐標
    az=np.zeros(n)           # az 標的物的方位角，北方為0，往東為正，正東方90度，正南方180度，正西方270度。
    alt=np.zeros(n)          # alt 標的物的仰角
    for i in range(n):
        coords[i]=np.array(np.vstack((e_e[:,i],e_n[:,i],e_z[:,i]))).dot(objs[i]-positions[i])
        r=np.sqrt(coords[i][0]**2+coords[i][1]**2)
        rho=np.sqrt(coords[i][0]**2+coords[i][1]**2+coords[i][2]**2)
        alt[i]=np.arcsin(coords[i][2]/rho)*180/np.pi
        az[i]=np.arctan2(coords[i][0],coords[i][1])*180/np.pi
        if az[i]<0:
            az[i]=az[i]+360
    if objname=="sun" or objname=="moon" :
      halfdisk=np.arcsin(r_obj/rho)*180/np.pi*u.deg
    else :
      halfdisk=0*u.deg
    thresh=-(halfdisk+Angle("0d34m0s")).deg   # 大氣折射，讓光線轉彎了34角分
    altdiff=alt-thresh        # altdif 變號的時候就是標的物升起或沉沒的區間，用內插進行估計
    for i in range(len(times)-1):
        if altdiff[i]*altdiff[i+1] <=0 :
            t=times[i]+(times[i+1]-times[i])*(0-altdiff[i])/(altdiff[i+1]-altdiff[i])
            if altdiff[i]<0:
                print("出",t,"方位",az[i])
            else:
                print("沒",t,"方位",az[i])
        if not (az[i] // 180 == az[i+1] // 180) and alt[i]>0:
            if np.abs(az[i]-180)<30:
                ratio=(180-az[i])/(az[i+1]-az[i])
                t=times[i]+(times[i+1]-times[i])*ratio
                print("中天",t,"仰角",alt[i]+(alt[i+1]-alt[i])*ratio,"方位","南")
            else:
                left=az[i]
                right=az[i+1]
                if left>300:
                    left=left-360
                if right>300:
                    right=right-360
                ratio=(0-left)/(az[i+1]-az[i])
                t=times[i]+(times[i+1]-times[i])*ratio
                print("中天",t,"仰角",alt[i]+(alt[i+1]-alt[i])*ratio,"方位","南")

if __name__ == "__main__":
    # execute only if run as a script
    import time
    time1 = time.time()
    # 用 multiprocessing 並不會比較快!
    # PrintCCalendar(2020)
    obj_sight0('2021-2-3','sun','0d32m0s','121.518d','25.063d',8)
    #obj_sight('2018-1-1','sun','121.5654d','25.032969d',8)
    time2 = time.time()
    print(time2-time1)


# 台北(東經121.5654度，北緯25.032969度，時差+8小時)，標的物太陽(sun)的視角32角分。
# Ex.1
# 用下面指令列印出2018年1月1日在台北的日出、日落、中天時刻以及方位與仰角。
# obj_sight('2018-1-1','sun','0d32m0s','121.5654d','25.032969d',8) 
# Rise 2018-01-01 06:38:47.130 Az 115.066472753
# Transit 2018-01-01 11:57:09.418 Alt 41.9575847767 Az South
# Set 2018-01-01 17:15:36.479 Az 244.828654564
# Ex.2
# obj_sight('2017-12-10','moon','0d32m0s','121.5654d','25.032969d',8)  
# Transit 2017-12-10 05:38:22.471 Alt 73.5230263955 Az South
# Set 2017-12-10 12:05:02.922 Az 278.440772592
# Ex.3
# obj_sight('2017-12-4','moon','0d32m0s','121.5654d','25.032969d',8)
# Set 2017-12-04 06:37:52.044 Az 290.09123571
# Rise 2017-12-04 17:54:55.105 Az 68.7458437964

