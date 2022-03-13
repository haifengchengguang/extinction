import numpy as np
import pandas as pd

df1 = pd.read_csv(r"allbdsdsswisegaia_notin_pos_fitsname_sdsswisegaia_extinc.csv")
rmag = df1['rmag']
imag = df1['imag']
zmag = df1['zmag']
Jmag = df1['Jmag']
Hmag = df1['Hmag']
Kmag = df1['Kmag']
W1mag = df1['W1mag']
W2mag = df1['W2mag']
W3mag = df1['W3mag']
W4mag = df1['W4mag']
extinc_mag = pd.concat(
        [rmag, imag, zmag, Jmag, Hmag, Kmag, W1mag, W2mag, W3mag, W4mag,rmag - imag, imag - zmag, zmag - Jmag, Jmag - Hmag,
         Hmag - Kmag, Kmag - W1mag, W1mag - W2mag, W2mag - W3mag, W3mag - W4mag, rmag - zmag, imag - Jmag, zmag - Hmag,
         Jmag - Kmag, Hmag - W1mag, Kmag - W2mag, W1mag - W3mag, W2mag - W4mag, rmag - Jmag, imag - Hmag, zmag - Kmag,
         Jmag - W1mag, Hmag - W2mag, Kmag - W3mag, W1mag - W4mag, rmag - Hmag, imag - Kmag, zmag - W1mag, Jmag - W2mag,
         Hmag - W3mag, Kmag - W4mag, rmag - Kmag, imag - W1mag, zmag - W2mag, Jmag - W3mag, Hmag - W4mag, rmag - W1mag,
         imag - W2mag, zmag - W3mag, Jmag - W4mag, rmag - W2mag, imag - W3mag, zmag - W4mag, rmag - W3mag, imag - W4mag,
         rmag - W4mag], axis=1)
np.savetxt('allbdsdsswisegaia_notin_pos_fitsname_sdsswisegaia_extinc_45.csv',extinc_mag,fmt="%f",delimiter=',')
print("end")