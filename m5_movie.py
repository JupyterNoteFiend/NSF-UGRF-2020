# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 18:50:39 2020

@author: Carleano Libretto
"""

from sunpy.net import Fido, attrs as a


from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = [10, 9]  # make plots larger

import sunpy.map
from sunpy.instr.aia import aiaprep
from sunpy.net import Fido, attrs as a
from astropy.coordinates import SkyCoord
from astropy import units as u
from matplotlib import animation
from IPython.display import HTML
import warnings
warnings.filterwarnings("ignore")

result = Fido.search(a.Time('2016/7/23 01:30:00', '2016/7/23 03:30:00'), a.Instrument.aia, a.Wavelength(304*u.Angstrom))[0,::25] 
#download = Fido.fetch(result, path ='C:\\Users\\Carleano Libretto\\downloads\\LIP2020\\{file}')
download = 'C:\\Users\\Carleano Libretto\\downloads\\LIP2020M5'
mapped_files = sunpy.map.Map(download)
aia_seq = []
k=0
for img in mapped_files:
    aiaprep(mapped_files[k])
    k+=1
    top_right = SkyCoord(1200*u.arcsec, 300*u.arcsec, frame=img.coordinate_frame)
    bottom_left = SkyCoord(600 * u.arcsec, -300. * u.arcsec, frame=img.coordinate_frame)
    aia_seq.append(img.submap(top_right, bottom_left))
aia = sunpy.map.Map(aia_seq)
fig, ax = plt.subplots()
# image half down sequence to get better scaling
plot_obj = aia[len(aia) // 2].plot()

def animate(i):
    ax.set_title("AIA %s %s" % (aia[i].meta['wave_str'],
                                aia[i].meta['t_obs'][:-8]))
    plot_obj.set_data(aia[i].data)
    return (plot_obj,)

anim = animation.FuncAnimation(fig, animate, init_func=None,
                               frames=len(aia), interval=500, blit=True)
anim
plt.close(fig)
HTML(anim.to_html5_video())
plt.show()