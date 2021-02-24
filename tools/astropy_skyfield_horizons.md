# Astropy，Skyfield，JPL Horizons 取得太陽系物體資料

* 基本上，[Astropy](https://docs.astropy.org/en/stable/) 的 get_obj 會得到從地心看obj的 GCRS apparent 赤道座標(gcrs_apparent)，對應於 [Skyfield](https://rhodesmill.org/skyfield/) 的 earth.at(t).observe(obj).apparent()。

* Astropy 用 gcrs_apparent.transform_to(GeocentricTrueEcliptic(equinox = t)) 將 GCRS apparent赤道座標轉成真時黃道座標(true ecliptic)，
對應於 Skyfield 的 gcrs_apparent.frame_latlon(ecliptic_frame)。

* Astropy 用 gcrs_apparent.transform_to(TETE(obstime = t)) 將 GCRS apparent赤道座標轉成真時赤道座標(True Ecliptic True Equator)，
對應於 Skyfield 的 gcrs_apparent.radec(epoch='date')。

* Astropy 用 gcrs_apparent.transform_to(AltAz(obstime=t, location=location))
 |將 GCRS apparent赤道座標轉成 AltAz　座標，對應於 Skyfield 的 gcrs_apparent.altaz()。註:Skyfield在get_body時已將觀測地資料存入。

* Skyfield 有 earth.at(t).observe(obj) 取得 GCRS astrometric 座標，Astropy 則沒有此功能。

* [JPL Horizons](https://ssd.jpl.nasa.gov/horizons.cgi) 則應有盡有最為齊全，但必須上網批次取得較為不便。JPL Horizons 可用 [astroquery.jplhorizons](https://astroquery.readthedocs.io/en/latest/jplhorizons/jplhorizons.html) 套件來程式檢索。

* Astrometric 只考慮了光線旅行的時間，Apparent 增加了光線重力彎曲與光線視差等因素。
* TETE,AzAlt,GeocentricTrueEcliptic 屬於動態座標系統。
  GCRS, BCRS 的坐標軸方向不隨時間變動，接近於TEME座標在J2000.0時的方向。
  