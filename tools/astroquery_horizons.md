# Horizons.ephemerides vs Horizons.vectors

Horizons.ephemerides 可以查得 RA, DEC, RA_app, DEC_app. Horizons.vectors 可以查得 x, y, z,
兩者之間有何關聯？

* Horizons.ephemerides 的時間尺度為 utc,  Horizons.vectors 則為 tdb
* Horizons.vectors(refplane='earth', aberrations='astrometric') 查得的x,y,z轉成球體座標就是RA/DEC
* Horizons.vectors(refplane='earth', aberrations='apparent') 查得的x,y,z，轉成TETE坐標系的球體座標卻與RA_app/DEC_app **有點差距**。
* Horizons.vectors(refplane='earth', aberrations='geometric') 所得的資料，經過光線時間的調整，可以轉成球體座標 RA/DEC

