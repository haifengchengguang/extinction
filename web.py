import numpy as np
import pandas as pd

import astropy.units as u
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarWebQuery
def distance_parallax(x):
    return abs(1000/x)
df1 = pd.read_csv(r"full_match_rizjhkw1_id_ra_dec_distance.csv")
df1head=df1.head(5)
ra_csv=df1['ra']
dec_csv=df1['dec_']
rmag_observe=df1['rmag']
imag_observe=df1['imag']
zmag_observe=df1['zmag']
Jmag_observe=df1['Jmag']
Hmag_observe=df1['Hmag']
Kmag_observe=df1['Kmag']
W1mag_observe=df1['W1mag']
W2mag_observe=df1['W2mag']
distance=df1['parallax'].map(distance_parallax)
ra_deg=ra_csv*u.deg
dec_deg=dec_csv*u.deg
distance_pc=distance*u.pc
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
coords = SkyCoord(ra_deg, dec_deg,distance=distance_pc, frame='icrs')

q = BayestarWebQuery(version='bayestar2015')
E = q(coords, mode='mean')
#imag=imag_o E*1.971
#print(E)
Epandas=pd.Series(E.tolist())
Epandas.fillna(0,inplace=True)
rmag=rmag_observe-Epandas*2.31
imag=imag_observe-Epandas*1.71
zmag=zmag_observe-Epandas*1.29
Jmag=Jmag_observe-Epandas*0.72
Hmag=Hmag_observe-Epandas*0.46
Kmag=Kmag_observe-Epandas*0.306
W1mag=W1mag_observe-Epandas*0.18
W2mag=W2mag_observe-Epandas*0.16
extinc_mag=pd.concat([rmag,imag,zmag,Jmag,Hmag,Kmag,W1mag,W2mag],axis=1)
#imag_ndarray=np.array(imag)
# print(type(E))
#print(E)
#print(imag_ndarray)
#np.savetxt('Extinction.csv',E,fmt="%f",delimiter=',')
#np.savetxt('imag_in.csv',imag_ndarray,fmt="%f",delimiter=',')
np.savetxt('total_mag_in.csv',extinc_mag,fmt="%f",delimiter=',')