#！！！这个代码只在linux上跑过，windows不行，mac没测试过
#pip install dustmaps
#网址
#https://dustmaps.readthedocs.io/en/latest/maps.html#bayestar
import numpy as np
import pandas as pd

import astropy.units as u
from astropy.coordinates import SkyCoord
from dustmaps.bayestar import BayestarWebQuery
def distance_parallax(x):
    return abs(1000/x)
df1 = pd.read_csv(r"lty_notin.csv")#这里改文件路径
#df1head=df1.head(5)
ra_csv=df1['ra']
dec_csv=df1['dec']

distance=df1['parallax'].map(distance_parallax)
ra_deg=ra_csv*u.deg
dec_deg=dec_csv*u.deg
distance_pc=distance*u.pc
#print(ra_deg)
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

q = BayestarWebQuery(version='bayestar2019')#这里要用baystar2019，2018应该也行，不能用2015,那个时候gaia数据还没出来，没适配
E = q(coords, mode='mean')
#imag=imag_o E*1.971
#print(E)
Epandas=pd.Series(E.tolist())
Epandas_bp_rp= 0.94 * ((0.901 * Epandas - 0.03) / 0.269) + 0.11#这个公式在这里https://iopscience.iop.org/article/10.3847/1538-4357/ab5362#apjab5362s4，4.5节 公式(30)(31)推导出来
Epandas_bp_rp.fillna(0, inplace=True)#Epandas_bp_rp是计算出来的消光，null值我置为了0，用gaia观测的bp-rp减去这个值即可
