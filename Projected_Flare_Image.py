# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 12:25:02 2020

@author: Carleano Libretto
"""
import warnings
warnings.simplefilter('ignore')
import sunpy
import numpy as np
from scipy.ndimage import gaussian_filter
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.colors
import astropy.units as u
import astropy.time
import astropy.constants as const
from astropy.coordinates import SkyCoord
from sunpy.map import GenericMap
from sunpy.util.metadata import MetaDict
from sunpy.sun import constants as sun_const
from sunpy.coordinates import Helioprojective,HeliographicStonyhurst,Heliocentric
from sunpy.instr.aia import aiaprep
import glob

def semi_circular_loop(length,theta0=0*u.deg):
    r_1 = const.R_sun
    def r_2_func(x):
        return np.arccos(0.5*x/r_1.to(u.cm).value) - np.pi + length.to(u.cm).value/2./x
    r_2 = scipy.optimize.bisect(r_2_func,length.to(u.cm).value/(2*np.pi),
                               length.to(u.cm).value/np.pi) * u.cm
    alpha = np.arccos(0.5*(r_2/r_1).decompose())
    phi = np.linspace(-np.pi*u.rad + alpha,np.pi*u.rad-alpha,2000)
    # Quadratic formula to find r
    a = 1.
    b = -2*(r_1.to(u.cm)*np.cos(phi.to(u.radian)))
    c = r_1.to(u.cm)**2 - r_2.to(u.cm)**2
    r = (-b + np.sqrt(b**2 - 4*a*c))/2/a 
    # Choose only points above the surface
    i_r = np.where(r>r_1)
    r = r[i_r]
    phi = phi[i_r]
    hcc_frame = Heliocentric(observer=SkyCoord(
        lon=75*u.deg,lat=25*u.deg,radius=r_1,frame='heliographic_stonyhurst'))
    return (SkyCoord(y=r.to(u.cm)*np.sin(phi.to(u.radian)),
                     x=u.Quantity(r.shape[0]*[0*u.cm]),
                     z=r.to(u.cm)*np.cos(phi.to(u.radian)),
                     frame=hcc_frame)
            .transform_to('heliographic_stonyhurst'))
loop = semi_circular_loop(140*u.Mm,theta0=5*u.deg)
aia = glob.glob('C:\\Users\\Carleano Libretto\\downloads\\LIP2020M5(171)')
aia = sorted(aia)
dummy_map = (sunpy.map.Map(aia))
dummy_map = aiaprep(dummy_map[0])
fig = plt.figure(figsize=(30,35))
ax = fig.gca(projection=dummy_map)
dummy_map.plot(title=False)
ax.plot_coord(loop.transform_to(dummy_map.coordinate_frame),color='C3',lw=3)
#ax.plot_coord(SkyCoord(1400*u.arcsec, 500*u.arcsec, frame=dummy_map.coordinate_frame))
#ax.plot_coord(SkyCoord(600*u.arcsec, -250*u.arcsec, frame=dummy_map.coordinate_frame))
