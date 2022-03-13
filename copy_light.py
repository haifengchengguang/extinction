# 首先需要安装三个包：
# astropy（pip install astropy），
# mwdust（pip install https://github.com/jobovy/mwdust）,
# isochrones（pip install isochrones）
# pip install mwdust==1.1
# 然后的步骤：
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord

import mwdust
from isochrones.mist.bc import MISTBolometricCorrectionGrid

bc_grid = MISTBolometricCorrectionGrid(["2MASS_J",
                                        "2MASS_H",
                                        "2MASS_Ks"])  # 热改正xxxx选MIST isochrones的热改正"Gaia_G_MAW"
# combined19A = mwdust.Combined19(filter='xxxx')  # 消光xxxx选取3d的G波段消光
# combined19Av = mwdust.Combined19(filter='CTIO V')
combined19AJ = mwdust.Combined19(filter='2MASS J')
combined19AH = mwdust.Combined19(filter='2MASS H')
combined19AKs = mwdust.Combined19(filter='2MASS Ks')
# sfd = mwdust.SFD(filter='')
Mbol_sun = 4.7554
Teff_sun = 5777


def get_Lbol(mag, distance, bc_band, A_band):
    return -0.4 * (mag - 5 * np.log10(distance / 10) + bc_band - A_band - Mbol_sun)


def get_Aks(ra, dec, distance):
    c_icrs = SkyCoord(ra=ra * u.degree, dec=dec * u.degree, frame='icrs')
    l = c_icrs.galactic.l.value
    b = c_icrs.galactic.b.value
    Aks = combined19AKs(l, b, distance / 1000)
    Aj = combined19AJ(l, b, distance / 1000)
    Ah = combined19AH(l, b, distance / 1000)
    return Aks[0], Aj[0], Ah[0]


def get_bc_ks(Teff, logg, feh, Av):
    return bc_grid.interp([Teff, logg, feh, Av])[0]


def get_L_bol(ra, dec, distance, Teff, logg, feh, Ks):
    A_ks, Av,Ah = get_Aks(ra, dec, distance)
    bc_ks = get_bc_ks(Teff, logg, feh, Av)
    log_Lbol2Lbol_sun = get_Lbol(Ks, distance, bc_ks, A_ks)

    return log_Lbol2Lbol_sun


if __name__ == '__main__':
    # Lbol=get_L_bol(ra, dec, distance, Teff, logg, feh, Ks)
    print('end')
