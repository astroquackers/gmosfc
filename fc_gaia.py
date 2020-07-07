from .fcutils import make_gaia
from .fcutils import writeregion
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


def gmos_longslit_gaiasyn(name, ra1, dec1, pm_ra_cosdec=0.0, pm_dec=0.0, pa=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75):


    """
    Function to create gmos finding charts from synthetic gaia images, given the target coordinates, proper motions and epoch of observation.
    Currently, the proper motion translation of the final finding charts is disabled, and the finding chart is produced in J2000 coordinate system.  However, the proper motion of the target is calculated for the observational epoch. If the new coordinates are significantly different to the (1st decimal degree) J2000 coordinates, the finding charts maybe inaccurate.
    
    Returns finding charts for single GMOS target at given coordinates.
    
    Parameters:
    name: String
    The name of the target.
    
    ra1: float
    Right Ascension of the target in J2000 Fk5 system.
    
    dec1: float
    Declination of the target in J2000 Fk5 system.
    
    pm_ra_cosdec: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_dec: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pa: float, optional
    Position angle of the target in decimal degrees. Must be between 0 and 360. Default is 0, not parallitic.
    
    slitwidth: float, optional
    Slit width of the GMOS slit in arcsec. Default is 1 arcsec.
    
    pnum: String, optional
    Gemini Program ID. Defaults to GMOS spectropscopy.
    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.


    
    """
##### Parsec input parameters
    time_obs = Time(time_obs, format='jyear_str')
    epoch = Time(epoch, format='jyear_str')
    
    pm_ra_cosdec = pm_ra_cosdec*u.mas/u.yr
    pm_dec = pm_dec = pm_dec*u.mas/u.yr
    coords = SkyCoord(ra1, dec1, unit=unit, pm_ra_cosdec=pm_ra_cosdec, pm_dec=pm_dec, frame=frame, obstime=epoch)
    coordsnew = coords.transform_to(FK5(equinox=time_obs))
    
    print ('The coordinates of', name, 'in the epoch', time_obs.value, 'are RA=',round(coordsnew.ra.value,5), ',DEC=',round(coordsnew.dec.value,5))
    
###### Make the GAIA image
###### note that currently, the epoch param is disabled, and set to default of J2000. Transforming coordinates from the GAIA DR2 catalogue without catalogue cleaning leads to weird images
###### Read noise is disabled. This can be set to the GMOS default values, however the QLT for the observer corrects for this noise, hence it is set to default. A poission noise image is instead added.
    gaia_image = make_gaia(coords, epoch=time_obs)

######
######
    print ('Creating finding chart using synthetic GAIA DR2 image for', name)
    

###### Create overlay of slit
###### Shift position angle by 90 to display in ds9.
###### pyregion is obselete
    pa = Angle(pa * u.deg)
    pa = pa.value
    slitwidth=slitwidth * u.arcsec
    slitwidths = (3.5*slitwidth.value)
    writeregion(name, coords.ra.value, coords.dec.value, slitwidths, pa+90)
    pa = round(pa,1)


######## Create figure
######## TeX instances removed for broaded compatability.
    fig = plt.figure()
    gc = aplpy.FITSFigure(gaia_image, north=True, figure=fig)
    gc.show_grayscale(stretch='linear', invert=True)
    gc.frame.set_linewidth(0.5)  # points
    gc.frame.set_color('black')
    gc.show_markers(coords.ra.value,coords.dec.value,edgecolor='red', facecolor="None", s=markersize, label='Target', lw=0.5)
    gc.show_regions(name + ".reg")
    #gc.show_rectangles([coords.ra.value],[coords.dec.value],0.09166666666666666,slitwidths, angle=pa, edgecolor='blue', lw=0.5, label='pa')
    gc.add_scalebar(1)
    gc.scalebar.show(0.03334, corner=corner, frame=True, borderpad=0.4, pad=0.5)
    gc.scalebar.set_corner(corner)
    gc.scalebar.set_frame(True)
    gc.scalebar.set_color('black')
    gc.scalebar.set_label("$2'$")
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_xtext('$\mathrm{Right\,Ascension (J2000)}$')
    gc.axis_labels.set_ytext('$\mathrm{Declination (J2000)}$')
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_font(size='medium', weight='medium', \
                        stretch='normal', family='serif', \
                        style='normal', variant='normal')
    gc.ticks.set_color("black")
    gc.tick_labels.set_xformat('hh:mm:ss')
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.tick_labels.set_font(size='small')
    gc.ticks.set_length(0)
    plt.tick_params(length=6, width=0.6, which='major', tickdir='in')
    plt.tick_params(length=3,  which='minor')
    gc.add_grid()
    gc.grid.show()
    gc.grid.set_alpha(0.5)
##### Generate compass radius.to(u.deg).value/np.cos(field_center.dec.radian)
    arrow_width = 20*u.arcsec
    radius = 165*u.arcsec
    x = coords.ra.value-radius.to(u.deg).value/np.cos(coords.dec.radian)
    y = coords.dec.value-radius.to(u.deg).value
    dx = arrow_width.to(u.deg).value
    dy = arrow_width.to(u.deg).value/np.cos(coords.dec.radian)
    
    gc.show_arrows(x,y,dx,0.0, color='black', width=2, head_length=10, head_width=10)
    gc.show_arrows(x,y,0.0,dy, color='black', width=2, head_length=10, head_width=10)

    fig.legend(ncol=1)
    gc.add_label(0.07, 0.95, "Gaia $G$", relative=True, color='black')
    title = str(pnum)+str('\n')+str(name)+' : Long-slit acquisition'
    fig.suptitle(title, fontsize=12, horizontalalignment='right')
    chartname = name + "_fc.jpg"
    gc.save(chartname, dpi=600, format="jpg")
    gc.close()
    plt.close()
    print ('Finding chart saved as',chartname)
    print ('Cleaning up')
    os.system("rm "+str(gaia_image))
    os.system("rm "+str(name)+".reg")



def gmos_blindoffset_gaiasyn(name, ra1, dec1, rablind, decblind, pm_ra_cosdec=0.0, pm_dec=0.0, pa=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75):


    """
    Function to create gmos finding charts from synthetic gaia images, given the target coordinates, proper motions and epoch of observation.
    Currently, the proper motion translation of the final finding charts is disabled, and the finding chart is produced in J2000 coordinate system.  However, the proper motion of the target is calculated for the observational epoch. If the new coordinates are significantly different to the (1st decimal degree) J2000 coordinates, the finding charts maybe inaccurate.
    
    Returns finding charts for blind offset GMOS target at given coordinates.
    
    Parameters:
    name: String
    The name of the target.
    
    ra1: float
    Right Ascension of the base target in J2000 Fk5 system.
    
    dec1: float
    Declination of the base target in J2000 Fk5 system.
    
    rablind: float
    Right Ascension of the blind offset target in J2000 Fk5 system.
    
    decblind: float
    Declination of the blind offset target in J2000 Fk5 system.

    
    pm_ra_cosdec: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_dec: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pa: float, optional
    Position angle of the target in decimal degrees. Must be between 0 and 360. Default is 0, not parallitic.
    
    slitwidth: float, optional
    Slit width of the GMOS slit in arcsec. Default is 1 arcsec.
    
    pnum: String, optional
    Gemini Program ID. Defaults to GMOS spectropscopy.
    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.


    
    """
##### Parsec input parameters
    time_obs = Time(time_obs, format='jyear_str')
    epoch = Time(epoch, format='jyear_str')
    
    pm_ra_cosdec = pm_ra_cosdec*u.mas/u.yr
    pm_dec = pm_dec = pm_dec*u.mas/u.yr
    coords = SkyCoord(ra1, dec1, unit=unit, pm_ra_cosdec=pm_ra_cosdec, pm_dec=pm_dec, frame=frame, obstime=epoch)
    coordsnew = coords.transform_to(FK5(equinox=time_obs))
    
    coordsblind = SkyCoord(rablind, decblind, unit=unit, frame=frame, obstime=epoch)
    coordsnewblind = coords.transform_to(FK5(equinox=time_obs))

##### Calculate separation
    
    sep = coordsblind.separation(coords)
    sep_arcsec = sep.arcsec
    if sep_arcsec > 300.:
        raise Exception ("A star too far! Please check your coordinates, as the separation between the Base and blind offset targets exceeds the recommended limit at {}".format(sep_arcsec))

    print ('The coordinates of', name, 'base target in the epoch', time_obs.value, 'are RA=',round(coordsnew.ra.value,5), ',DEC=',round(coordsnew.dec.value,5))
    
    print ('The coordinates of', name, 'in the epoch', time_obs.value, 'are RA=',round(coordsnewblind.ra.value,5), ',DEC=',round(coordsnewblind.dec.value,5))

    
###### Make the GAIA image
###### note that currently, the epoch param is disabled, and set to default of J2000. Transforming coordinates from the GAIA DR2 catalogue without catalogue cleaning leads to weird images
###### Read noise is disabled. This can be set to the GMOS default values, however the QLT for the observer corrects for this noise, hence it is set to default. A poission noise image is instead added.
    gaia_image = make_gaia(coords, epoch=time_obs)

######
######
    print ('Creating finding chart using synthetic GAIA DR2 image for', name)
    

###### Create overlay of slit
###### Shift position angle by 90 to display in ds9.
###### pyregion is obselete
    pa = Angle(pa * u.deg)
    pa = pa.value
    slitwidth=slitwidth * u.arcsec
    slitwidths = (3.5*slitwidth.value)
    writeregion(name, coords.ra.value, coords.dec.value, slitwidths, pa+90)
    pa = round(pa,1)


######## Create figure
######## TeX instances removed for broaded compatability.
    fig = plt.figure()
    gc = aplpy.FITSFigure(gaia_image, north=True, figure=fig)
    gc.show_grayscale(stretch='linear', invert=True)
    gc.frame.set_linewidth(0.5)  # points
    gc.frame.set_color('black')
    gc.show_markers(coords.ra.value,coords.dec.value,edgecolor='red', facecolor="None", s=markersize, label='Base', lw=0.5)
    gc.show_regions(name + ".reg")
    gc.show_markers(coordsblind.ra.value,coordsblind.dec.value,edgecolor='blue', facecolor="None", marker='s', s=markersize, label='Blind offset', lw=0.5)
 #gc.show_rectangles([coords.ra.value],[coords.dec.value],0.09166666666666666,slitwidths, angle=pa, edgecolor='blue', lw=0.5, label='pa')
    gc.add_scalebar(1)
    gc.scalebar.show(0.03334, corner=corner, frame=True, borderpad=0.4, pad=0.5)
    gc.scalebar.set_corner(corner)
    gc.scalebar.set_frame(True)
    gc.scalebar.set_color('black')
    gc.scalebar.set_label("$2'$")
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_xtext('$\mathrm{Right\,Ascension (J2000)}$')
    gc.axis_labels.set_ytext('$\mathrm{Declination (J2000)}$')
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_font(size='medium', weight='medium', \
                        stretch='normal', family='serif', \
                        style='normal', variant='normal')
    gc.ticks.set_color("black")
    gc.tick_labels.set_xformat('hh:mm:ss')
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.tick_labels.set_font(size='small')
    gc.ticks.set_length(0)
    plt.tick_params(length=6, width=0.6, which='major', tickdir='in')
    plt.tick_params(length=3,  which='minor')
    gc.add_grid()
    gc.grid.show()
    gc.grid.set_alpha(0.5)
##### Generate compass radius.to(u.deg).value/np.cos(field_center.dec.radian)
    arrow_width = 20*u.arcsec
    radius = 165*u.arcsec
    x = coords.ra.value-radius.to(u.deg).value/np.cos(coords.dec.radian)
    y = coords.dec.value-radius.to(u.deg).value
    dx = arrow_width.to(u.deg).value
    dy = arrow_width.to(u.deg).value/np.cos(coords.dec.radian)
    
    gc.show_arrows(x,y,dx,0.0, color='black', width=2, head_length=10, head_width=10)
    gc.show_arrows(x,y,0.0,dy, color='black', width=2, head_length=10, head_width=10)

    fig.legend(ncol=1)
    gc.add_label(0.07, 0.95, "Gaia $G$", relative=True, color='black')
    title = str(pnum)+str('\n')+str(name)+' : Blind-offset acquisition'
    fig.suptitle(title, fontsize=12, horizontalalignment='right')
    chartname = name + "_fc.jpg"
    gc.save(chartname, dpi=600, format="jpg")
    gc.close()
    plt.close()
    print ('Finding chart saved as',chartname)
    print ('Cleaning up')
    os.system("rm "+str(gaia_image))
    os.system("rm "+str(name)+".reg")
















def gmos_twotarget_gaiasyn(name, ra1, dec1, ra2, dec2, pm_ra_cosdec1=0.0, pm_dec1=0.0, pm_ra_cosdec2=0.0, pm_dec2=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75, markersizecentral=25):


    """
    Function to create gmos finding charts from synthetic gaia images, given the target coordinates, proper motions and epoch of observation.
    Currently, the proper motion translation of the final finding charts is disabled, and the finding chart is produced in J2000 coordinate system.  However, the proper motion of the target is calculated for the observational epoch. If the new coordinates are significantly different to the (1st decimal degree) J2000 coordinates, the finding charts maybe inaccurate.
    
    Returns finding charts for two GMOS targets at given coordinates.
    
    Parameters:
    name: String
    The name of the target.
    
    ra1: float
    Right Ascension of the user 1 in J2000 Fk5 system.
    
    dec1: float
    Declination of the user 1 in J2000 Fk5 system.
    
    ra2: float
    Right Ascension of the user 2 in J2000 Fk5 system.
    
    dec2: float
    Declination of the user 2 in J2000 Fk5 system.

    
    pm_ra_cosdec1: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_dec1: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_ra_cosdec2: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_dec2: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    slitwidth: float, optional
    Slit width of the GMOS slit in arcsec. Default is 1 arcsec.
    
    pnum: String, optional
    Gemini Program ID. Defaults to GMOS spectropscopy.
    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.


    
    """
##### Parsec input parameters
    time_obs = Time(time_obs, format='jyear_str')
    epoch = Time(epoch, format='jyear_str')
    
    pm_ra_cosdec1 = pm_ra_cosdec1*u.mas/u.yr
    pm_dec1 = pm_dec1*u.mas/u.yr
    coords1 = SkyCoord(ra1, dec1, unit=unit, pm_ra_cosdec=pm_ra_cosdec1, pm_dec=pm_dec1, frame=frame, obstime=epoch)
    coordsnew1 = coords1.transform_to(FK5(equinox=time_obs))
    
    pm_ra_cosdec2 = pm_ra_cosdec2*u.mas/u.yr
    pm_dec2 = pm_dec2*u.mas/u.yr
    coords2 = SkyCoord(ra2, dec2, unit=unit, pm_ra_cosdec=pm_ra_cosdec2, pm_dec=pm_dec2, frame=frame, obstime=epoch)
    coordsnew2 = coords2.transform_to(FK5(equinox=time_obs))

##### Calculate separation
    
    sep = coords2.separation(coords1)
    sep_arcsec = sep.arcsec
    if sep_arcsec > 330.:
        raise Exception ("A star too far! Please check your coordinates, as the separation between User 1 and User 2 targets exceeds the recommended limit at {}".format(sep_arcsec))

    print ('The coordinates of', name, 'User 1 target in the epoch', time_obs.value, 'are RA=',round(coordsnew1.ra.value,5), ',DEC=',round(coordsnew1.dec.value,5))
    
    print ('The coordinates of', name, 'User 2 target in the epoch', time_obs.value, 'are RA=',round(coordsnew2.ra.value,5), ',DEC=',round(coordsnew2.dec.value,5))


###### Calculate position angle and Base coordinates

    meanra = (coords1.ra.to_value()+coords2.ra.to_value())/2
    meandec = (coords1.dec.to_value()+coords2.dec.to_value())/2
    
    coords = SkyCoord(meanra, meandec, unit=unit, frame=frame, obstime=epoch)
    
    pa = coords1.position_angle(coords2).degree

    log.info("The position angle is %+.2f" % pa )
    log.info("The Base target in Right Ascencsion (J2000) is %+.2f" % meanra)
    log.info("The Base target in Declination is (J2000) is %+.2f" % meandec )

###### Make the GAIA image
###### note that currently, the epoch param is disabled, and set to default of J2000. Transforming coordinates from the GAIA DR2 catalogue without catalogue cleaning leads to weird images
###### Read noise is disabled. This can be set to the GMOS default values, however the QLT for the observer corrects for this noise, hence it is set to default. A poission noise image is instead added.
    gaia_image = make_gaia(coords, epoch=time_obs)

######
######
    print ('Creating finding chart using synthetic GAIA DR2 image for', name)
    

###### Create overlay of slit
###### Shift position angle by 90 to display in ds9.
###### pyregion is obselete
    pa = Angle(pa * u.deg)
    pa = pa.value
    slitwidth=slitwidth * u.arcsec
    slitwidths = (3.5*slitwidth.value)
    writeregion(name, coords.ra.value, coords.dec.value, slitwidths, pa+90)
    pa = round(pa,1)


######## Create figure
######## TeX instances removed for broaded compatability.
    fig = plt.figure()
    gc = aplpy.FITSFigure(gaia_image, north=True, figure=fig)
    gc.show_grayscale(stretch='linear', invert=True)
    gc.frame.set_linewidth(0.5)  # points
    gc.frame.set_color('black')
    gc.show_markers(coords.ra.value,coords.dec.value,edgecolor='red', facecolor="None", marker='o', s=markersizecentral, label='Base', lw=0.5)

    gc.show_markers(coords1.ra.value,coords1.dec.value,edgecolor='blue', facecolor="None", marker='s', s=markersize, label='User1', lw=0.5)
    gc.show_regions(name + ".reg")
    gc.show_markers(coords2.ra.value,coords2.dec.value,edgecolor='green', facecolor="None", marker='s', s=markersize, label='User2', lw=0.5)
 #gc.show_rectangles([coords.ra.value],[coords.dec.value],0.09166666666666666,slitwidths, angle=pa, edgecolor='blue', lw=0.5, label='pa')
    gc.add_scalebar(1)
    gc.scalebar.show(0.03334, corner=corner, frame=True, borderpad=0.4, pad=0.5)
    gc.scalebar.set_corner(corner)
    gc.scalebar.set_frame(True)
    gc.scalebar.set_color('black')
    gc.scalebar.set_label("$2'$")
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_xtext('$\mathrm{Right\,Ascension (J2000)}$')
    gc.axis_labels.set_ytext('$\mathrm{Declination (J2000)}$')
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_font(size='medium', weight='medium', \
                        stretch='normal', family='serif', \
                        style='normal', variant='normal')
    gc.ticks.set_color("black")
    gc.tick_labels.set_xformat('hh:mm:ss')
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.tick_labels.set_font(size='small')
    gc.ticks.set_length(0)
    plt.tick_params(length=6, width=0.6, which='major', tickdir='in')
    plt.tick_params(length=3,  which='minor')
    gc.add_grid()
    gc.grid.show()
    gc.grid.set_alpha(0.5)
##### Generate compass radius.to(u.deg).value/np.cos(field_center.dec.radian)
    arrow_width = 20*u.arcsec
    radius = 165*u.arcsec
    x = coords.ra.value-radius.to(u.deg).value/np.cos(coords.dec.radian)
    y = coords.dec.value-radius.to(u.deg).value
    dx = arrow_width.to(u.deg).value
    dy = arrow_width.to(u.deg).value/np.cos(coords.dec.radian)
    
    gc.show_arrows(x,y,dx,0.0, color='black', width=2, head_length=10, head_width=10)
    gc.show_arrows(x,y,0.0,dy, color='black', width=2, head_length=10, head_width=10)

    fig.legend(ncol=1)
    gc.add_label(0.07, 0.95, "Gaia $G$", relative=True, color='black')
    title = str(pnum)+str('\n')+str(name)+' : Two target acquisition'
    fig.suptitle(title, fontsize=12, horizontalalignment='right')
    chartname = name + "_fc.jpg"
    gc.save(chartname, dpi=600, format="jpg")
    gc.close()
    plt.close()
    print ('Finding chart saved as',chartname)
    print ('Cleaning up')
    os.system("rm "+str(gaia_image))
    os.system("rm "+str(name)+".reg")





#####################
#####################
########## SKY SURVEY FINDING CHARTS



def gmos_twotarget(name, ra1, dec1, ra2, dec2, pm_ra_cosdec1=0.0, pm_dec1=0.0, pm_ra_cosdec2=0.0, pm_dec2=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', skysurvey="DSS", frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75, markersizecentral=25):


    """
    Function to create gmos finding charts from archive sky images, given the target coordinates, proper motions and epoch of observation.
    Currently, the proper motion translation of the final finding charts is disabled, and the finding chart is produced in J2000 coordinate system.  However, the proper motion of the target is calculated for the observational epoch. If the new coordinates are significantly different to the (1st decimal degree) J2000 coordinates, the finding charts maybe inaccurate.
    
    Returns finding charts for two GMOS targets at given coordinates.
    
    Parameters:
    name: String
    The name of the target.
    
    ra1: float
    Right Ascension of the user 1 in J2000 Fk5 system.
    
    dec1: float
    Declination of the user 1 in J2000 Fk5 system.
    
    ra2: float
    Right Ascension of the user 2 in J2000 Fk5 system.
    
    dec2: float
    Declination of the user 2 in J2000 Fk5 system.

    
    pm_ra_cosdec1: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_dec1: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_ra_cosdec2: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    pm_dec2: float, optional
    Proper motion of the target in mas/yr. Depreciated.
    
    slitwidth: float, optional
    Slit width of the GMOS slit in arcsec. Default is 1 arcsec.
    
    pnum: String, optional
    Gemini Program ID. Defaults to GMOS spectropscopy.
    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.


    
    """
##### Parsec input parameters
    time_obs = Time(time_obs, format='jyear_str')
    epoch = Time(epoch, format='jyear_str')
    
    pm_ra_cosdec1 = pm_ra_cosdec1*u.mas/u.yr
    pm_dec1 = pm_dec1*u.mas/u.yr
    coords1 = SkyCoord(ra1, dec1, unit=unit, pm_ra_cosdec=pm_ra_cosdec1, pm_dec=pm_dec1, frame=frame, obstime=epoch)
    coordsnew1 = coords1.transform_to(FK5(equinox=time_obs))
    
    pm_ra_cosdec2 = pm_ra_cosdec2*u.mas/u.yr
    pm_dec2 = pm_dec2*u.mas/u.yr
    coords2 = SkyCoord(ra2, dec2, unit=unit, pm_ra_cosdec=pm_ra_cosdec2, pm_dec=pm_dec2, frame=frame, obstime=epoch)
    coordsnew2 = coords2.transform_to(FK5(equinox=time_obs))

##### Calculate separation
    
    sep = coords2.separation(coords1)
    sep_arcsec = sep.arcsec
    if sep_arcsec > 330.:
        raise Exception ("A star too far! Please check your coordinates, as the separation between User 1 and User 2 targets exceeds the recommended limit at {}".format(sep_arcsec))

    print ('The coordinates of', name, 'User 1 target in the epoch', time_obs.value, 'are RA=',round(coordsnew1.ra.value,5), ',DEC=',round(coordsnew1.dec.value,5))
    
    print ('The coordinates of', name, 'User 2 target in the epoch', time_obs.value, 'are RA=',round(coordsnew2.ra.value,5), ',DEC=',round(coordsnew2.dec.value,5))


###### Calculate position angle and Base coordinates

    meanra = (coords1.ra.to_value()+coords2.ra.to_value())/2
    meandec = (coords1.dec.to_value()+coords2.dec.to_value())/2
    
    coords = SkyCoord(meanra, meandec, unit=unit, frame=frame, obstime=epoch)
    
    pa = coords1.position_angle(coords2).degree

    log.info("The position angle is %+.2f" % pa )
    log.info("The Base target in Right Ascencsion (J2000) is %+.2f" % meanra)
    log.info("The Base target in Declination is (J2000) is %+.2f" % meandec )

###### Create the sky survey image
    imageurl = SkyView.get_image_list(position=coords, survey=skysurvey, coordinates="J2000", height=7*u.arcmin, width=7*u.arcmin)
    filename = wget.download(imageurl[0])
    gaia_image = str(filename)

######
######
    print ('Creating finding chart using archival survey imaging for', name)
    

###### Create overlay of slit
###### Shift position angle by 90 to display in ds9.
###### pyregion is obselete
    pa = Angle(pa * u.deg)
    pa = pa.value
    slitwidth=slitwidth * u.arcsec
    slitwidths = (3.5*slitwidth.value)
    writeregion(name, coords.ra.value, coords.dec.value, slitwidths, pa+90)
    pa = round(pa,1)


######## Create figure
######## TeX instances removed for broaded compatability.
    fig = plt.figure()
    gc = aplpy.FITSFigure(gaia_image, north=True, figure=fig)
    gc.show_grayscale(stretch='linear', invert=True)
    gc.frame.set_linewidth(0.5)  # points
    gc.frame.set_color('black')
    gc.show_markers(coords.ra.value,coords.dec.value,edgecolor='red', facecolor="None", marker='o', s=markersizecentral, label='Base', lw=0.5)

    gc.show_markers(coords1.ra.value,coords1.dec.value,edgecolor='blue', facecolor="None", marker='s', s=markersize, label='User1', lw=0.5)
    gc.show_regions(name + ".reg")
    gc.show_markers(coords2.ra.value,coords2.dec.value,edgecolor='green', facecolor="None", marker='s', s=markersize, label='User2', lw=0.5)
 #gc.show_rectangles([coords.ra.value],[coords.dec.value],0.09166666666666666,slitwidths, angle=pa, edgecolor='blue', lw=0.5, label='pa')
    gc.add_scalebar(1)
    gc.scalebar.show(0.03334, corner=corner, frame=True, borderpad=0.4, pad=0.5)
    gc.scalebar.set_corner(corner)
    gc.scalebar.set_frame(True)
    gc.scalebar.set_color('black')
    gc.scalebar.set_label("$2'$")
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_xtext('$\mathrm{Right\,Ascension (J2000)}$')
    gc.axis_labels.set_ytext('$\mathrm{Declination (J2000)}$')
    gc.scalebar.set_font(size='small', weight='medium', \
                     stretch='normal', family='serif', \
                     style='normal', variant='normal')
    gc.axis_labels.set_font(size='medium', weight='medium', \
                        stretch='normal', family='serif', \
                        style='normal', variant='normal')
    gc.ticks.set_color("black")
    gc.tick_labels.set_xformat('hh:mm:ss')
    gc.tick_labels.set_yformat('dd:mm:ss')
    gc.tick_labels.set_font(size='small')
    gc.ticks.set_length(0)
    plt.tick_params(length=6, width=0.6, which='major', tickdir='in')
    plt.tick_params(length=3,  which='minor')
    gc.add_grid()
    gc.grid.show()
    gc.grid.set_alpha(0.5)
##### Generate compass radius.to(u.deg).value/np.cos(field_center.dec.radian)
    arrow_width = 20*u.arcsec
    radius = 165*u.arcsec
    x = coords.ra.value-radius.to(u.deg).value/np.cos(coords.dec.radian)
    y = coords.dec.value-radius.to(u.deg).value
    dx = arrow_width.to(u.deg).value
    dy = arrow_width.to(u.deg).value/np.cos(coords.dec.radian)
    
    gc.show_arrows(x,y,dx,0.0, color='black', width=2, head_length=5, head_width=5)
    gc.show_arrows(x,y,0.0,dy, color='black', width=2, head_length=5, head_width=5)

    fig.legend(ncol=1)
    gc.add_label(0.07, 0.95, skysurvey, relative=True, color='black')
    title = str(pnum)+str('\n')+str(name)+' : Two target acquisition'
    fig.suptitle(title, fontsize=12, horizontalalignment='right')
    chartname = name + "_fc.jpg"
    gc.save(chartname, dpi=600, format="jpg")
    gc.close()
    plt.close()
    print ('Finding chart saved as',chartname)
    print ('Cleaning up')
    os.system("rm "+str(gaia_image))
    os.system("rm "+str(name)+".reg")





