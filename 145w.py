import numpy as np
import pandas as pd

import astropy.units as u
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarWebQuery
def distance_parallax(x):
    return abs(1000/x)
df1 = pd.read_csv(r"lowmass_sdss_2mass_wise_gaia.csv",nrows=800000)
#df_2=pd.read_csv(r"lowmass_sdss_2mass_wise_gaia.csv",skiprows=800000)


print("read file")
ra_csv = df1['ra']
dec_csv = df1['dec_']
rmag_observe = df1['rmag']
imag_observe = df1['imag']
zmag_observe = df1['zmag']
Jmag_observe = df1['Jmag']
Hmag_observe = df1['Hmag']
Kmag_observe = df1['Kmag']
W1mag_observe = df1['W1mag']
W2mag_observe = df1['W2mag']
W3mag = df1['W3mag']
W4mag = df1['W4mag']
e_rmag = df1['e_rmag']
e_imag = df1['e_imag']
e_zmag = df1['e_zmag']
e_Jmag = df1['e_Jmag']
e_Hmag = df1['e_Hmag']
e_Kmag = df1['e_Kmag']
e_W1mag = df1['e_W1mag']
e_W2mag = df1['e_W2mag']
e_W3mag = df1['e_W3mag']
e_W4mag = df1['e_W4mag']
id = df1['id']
ab_i = df1['ad_i']
parallax = df1['parallax']
distance_ = df1['distance']
pm = df1['pm']
distance = df1['parallax'].map(distance_parallax)
ra_deg = ra_csv * u.deg
dec_deg = dec_csv * u.deg
distance_pc = distance * u.pc
# b0 = [10., 12., -25.]
# l0=[90., 150., 35.]
# data={
#     "ra":l0,
#     "dec":b0,
# }
# df=pd.DataFrame(data)
# myvar_l=df['ra']
# l = myvar_l* u.deg
# myvar_b=df['dec']
# b = myvar_b* u.deg
# d = [500., 3500., 1000.] * u.pc
# print(type(l0))
# distance=d,
print("variable")
coords = SkyCoord(ra_deg, dec_deg, distance=distance_pc, frame='icrs')
q = BayestarWebQuery(version='bayestar2015')
E = q(coords, mode='mean')
# imag=imag_o E*1.971
# print(E)
print("dustmap")
Epandas = pd.Series(E.tolist())
Epandas.fillna(0, inplace=True)
rmag = rmag_observe - Epandas * 2.31
imag = imag_observe - Epandas * 1.71
zmag = zmag_observe - Epandas * 1.29
Jmag = Jmag_observe - Epandas * 0.72
Hmag = Hmag_observe - Epandas * 0.46
Kmag = Kmag_observe - Epandas * 0.306
W1mag = W1mag_observe - Epandas * 0.18
W2mag = W2mag_observe - Epandas * 0.16
print("calculate")
extinc_mag = pd.concat(
        [rmag, imag, zmag, Jmag, Hmag, Kmag, W1mag, W2mag, W3mag, W4mag], axis=1)
    # imag_ndarray=np.array(imag)

# ra_csv, dec_csv, , e_rmag, e_imag, e_zmag,
#          e_Jmag, e_Hmag, e_Kmag, e_W1mag, e_W2mag, e_W3mag, e_W4mag, rmag - imag, imag - zmag, zmag - Jmag, Jmag - Hmag,
#          Hmag - Kmag, Kmag - W1mag, W1mag - W2mag, W2mag - W3mag, W3mag - W4mag, rmag - zmag, imag - Jmag, zmag - Hmag,
#          Jmag - Kmag, Hmag - W1mag, Kmag - W2mag, W1mag - W3mag, W2mag - W4mag, rmag - Jmag, imag - Hmag, zmag - Kmag,
#          Jmag - W1mag, Hmag - W2mag, Kmag - W3mag, W1mag - W4mag, rmag - Hmag, imag - Kmag, zmag - W1mag, Jmag - W2mag,
#          Hmag - W3mag, Kmag - W4mag, rmag - Kmag, imag - W1mag, zmag - W2mag, Jmag - W3mag, Hmag - W4mag, rmag - W1mag,
#          imag - W2mag, zmag - W3mag, Jmag - W4mag, rmag - W2mag, imag - W3mag, zmag - W4mag, rmag - W3mag, imag - W4mag,
#          rmag - W4mag, id, ab_i, parallax, distance_, pm

#df1head=df1.head(5)
# extinc_mag_1=dustmapsql(df_1)
# extinc_mag_2=dustmapsql(df_2)
# extinc_mag=pd.concat([extinc_mag_1,extinc_mag_2],axis=1)
# print(type(E))
#print(E)
#print(imag_ndarray)
#np.savetxt('Extinction.csv',E,fmt="%f",delimiter=',')
#np.savetxt('imag_in.csv',imag_ndarray,fmt="%f",delimiter=',')
np.savetxt('lowmass_sdss_2mass_wise_gaia_extinc.csv',extinc_mag,fmt="%f",delimiter=',')