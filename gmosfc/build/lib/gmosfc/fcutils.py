import matplotlib.pyplot as plt
import numpy as np
import os
from astroquery.skyview import SkyView
import wget
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5
from astropy.coordinates import Angle, Latitude, Longitude
import astropy.units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.time import Time
from astropy import log
from astropy.wcs import WCS
from astroquery.gaia import Gaia
from astropy.table import QTable, Table, Column
from astropy.time import Time
from photutils.datasets import make_gaussian_sources_image
from photutils.datasets import make_noise_image

import aplpy







def writeregion(target, ra,dec,slitwidth,pa):
    file1 = open(target + ".reg","w")
    paround = round((pa-90),1)
    newstr = str('box(')+str(ra)+str(',')+str(dec)+str(',330.000",')+str(slitwidth)+str('",')+str(pa)+str(')'+str(' # text={')+str(paround)+str('} dash=1'))
    L = ['global color=blue dashlist=8 3 width=1 font="helvetica 7 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n','fk5 \n', newstr]
    file1.writelines(L)
    file1.close()



    


    




def make_gaia(coords, image=os.path.join('.','synthetic_gaia.fits'), epoch="J2020.5"):
    """
    Creates a synthetic image from GAIA DR2 photometry along the specified coordinates.

    Returns the synthetic image in fits format.

    Parameters:
    coords: Coordinates of centre of image.
    image: File name of image. Default is synthetic_gaia.fits
    epoch: Epoch to translate image. Default is J2020.5
    """
#####
##### Define GMOS parameters
##### GMOS FoV is 330 arcsec
##### GMOS slit length is 330 arcsec
##### GMOS fwhm chosen in 0.3 arcsec
##### Image created is 390 arcsec, to accomadte additional overlays
##### GMOS read noise is 3.96  (e-/ADU)
##### GMOS gain is 1.829  (e- rms)
#####

    gmosfov = 390 * u.arcsec
    fwhm = 0.3 * u.arcsec
    shape = (1300,1300)
    zeroarr = np.zeros(shape)
    width = 390 * u.arcsec
    height = 390 * u.arcsec
    
##### Read coordinates and epoch
    coords=coords
    epoch = Time(epoch, format='jyear_str')
    
##### Query GAIA DR2 using astroquery TAP+
##### Synchronous query, limit of 2000 rows
    print ('Searching GAIA TAP+ for stars over the GMOS FoV...')
    gaiasearch = Gaia.cone_search_async(coordinate=coords, radius=gmosfov, verbose=False)
    searchresults = gaiasearch.get_results()
    nstars = len(searchresults)
    print(nstars, 'stars recovered from GAIA TAP+ query over the GMOS FoV')
    if nstars == 2000:
       warnings.warn('Asynchronous TAP+ service query limit reached, image is incomplete.')


    fitshdr = mkfitshdr(coords, zeroarr, image, epoch)
    hdr = fitshdr[0]
    hdu = fitshdr[1]
    hdul = fitshdr[2]
    w = fitshdr[3]
    
    c = SkyCoord(ra=searchresults['ra'] ,dec=searchresults['dec'], pm_ra_cosdec=searchresults['pmra'], pm_dec=searchresults['pmdec'], obstime=Time(searchresults['ref_epoch'], format='decimalyear'))
    
    cnew = c
###### in the current implementation, transformation of the image to the epoch is not applied

######    cnew = c.apply_space_motion(new_obstime=Time(epoch))
    star_coords = w.wcs_world2pix([cnew.ra.deg],[cnew.dec.deg],0)
   
   
    xvalue = star_coords[0].ravel()
    yvalue = star_coords[1].ravel()
    flux = (searchresults['phot_g_mean_flux'] * u.mag).value

##### Error values
    cerr = SkyCoord(ra=searchresults['ra_error'] ,dec=searchresults['dec_error'])

    xerrvalue = cerr.ra.value
    yerrvalue = cerr.dec.value
    
    xerrmask = np.nan_to_num(xerrvalue, nan=0.0, posinf=0.0, neginf=0.0)
    yerrmask = np.nan_to_num(yerrvalue, nan=0.0, posinf=0.0, neginf=0.0)
    xerrmask[xerrmask < 1.0] = 1.0
    yerrmask[yerrmask < 1.0] = 1.0
##### Error values less than 1.0 lead to incorrect values in photutils. Error values of 1.0 lead to no significiant differences in images compared to no error.

##### Create Table of results


    t = Table([flux, xvalue, yvalue, xerrmask, yerrmask], names=('flux','x_mean','y_mean', 'x_stddev', 'y_stddev'), dtype=('i8', 'i8', 'i8', 'i8', 'i8'))
    print ('Table of stars is being generated as an image')
    zeroarr = make_gaussian_sources_image(shape,t)

##### Read noise

    read_noise = np.random.normal(scale=(3.92/1.829), size=shape)
    noise_image = make_noise_image(shape, distribution='poisson', mean=0.5)
#####

    synth_image = zeroarr + noise_image
#####
###### read_noise is disabled
    hdu = fits.PrimaryHDU(synth_image, header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.writeto(image, overwrite=True)
    print ('Synthetic GAIA DR2 image is saved as',image)
    return image




def mkfitshdr (coords, zeroarr, image, epoch):
    """
    Creates FITS header for synthetic GAIA images
    
    Returns FITS header and WCS
    
    Parameters:
    coords: Central coordinates of fits image
    zeroarr: Shape of image
    image: Image file name
    """
#### Define image fwhm
    fwhm = 0.3 * u.arcsec
    epoch = epoch.value
    ####### Curently epoch is disabled, and set by default to 2020
    hdr = fits.Header()
    hdr.set('CRPIX1', 650, 'X reference pixel')
    hdr.set('CRPIX2', 650, 'Y reference pixel')
    hdr.set('CRVAL1', coords.ra.deg, 'Reference longitude')
    hdr.set('CRVAL2', coords.dec.deg, 'Reference latitude')
    hdr.set('CTYPE1', 'RA---TAN', 'Coordinates -- projection')
    hdr.set('CTYPE2', 'DEC--TAN', 'Coordinates -- projection')
    hdr.set('CDELT1', -fwhm.to(u.deg).value, 'X scale')# mind the sign !
    hdr.set('CDELT2', fwhm.to(u.deg).value, 'Y scale')
    hdr.set('RADESYS', 'FK5', 'Coordinate system')
    hdr.set('EQUINOX', 'J2000', 'Epoch of the equinox')
    hdu = fits.PrimaryHDU(zeroarr, header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.writeto(image, overwrite=True)
    w = WCS(image)
    return hdr, hdu, hdul, w

