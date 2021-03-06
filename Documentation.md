Installation
----------
gmosfc can be installed via pip using the command

    python3 -m pip install gmosfc
    
gmosfc requires python3 to run. It will not work with the python2. Please also ensure the dependencies are installed before installing gmosfc. The current version (beta) is 0.0.11.
Installation has been tested on MAC and LINUX (opensuse) systems. If it has a glitch on your system, please let us know.

Reference guide
----------

gmosfc includes functinos to create finding charts either from the astroquery available imaging sky surveys, or using synthetic GAIA DR2 images. There are six functions, each are described here.

Before beginning, we need to import gmosfc

    import gmosfc as gfc
    
We can now run functions avaiable in gmosfc as
 
    gfc.gmos_blindoffset()
    
All coordinates used are in astropy SkyCoord format, described here https://docs.astropy.org/en/stable/api/astropy.coordinates.SkyCoord.html
    
The available functions for gfc using synthetic GAIA DR2 images are

1. Longslit GMOS target with synthetic GAIA images

```python
gmos_longslit_gaiasyn(name, ra1, dec1, pm_ra_cosdec=0.0, pm_dec=0.0, pa=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75):

    """
    Function to create gmos finding charts from synthetic gaia images, given the target coordinates, proper motions and epoch of observation.
    Currently, the proper motion translation of the final finding charts is disabled, and the finding chart is produced in J2000 coordinate system.  However,  the proper motion of the target is calculated for the observational epoch. If the new coordinates are significantly different to the (1st decimal degree) J2000 coordinates, the finding charts maybe inaccurate.
    
    Returns finding charts for single GMOS target at given coordinates.
    
    Parameters:
    name: String
    The name of the target.
    
    ra1: float
    Right Ascension of the target in J2000 FK5 system. In astropy system
    
    dec1: float
    Declination of the target in J2000 FK5 system. In astropy system.
    
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

    markersize: float, optional
    Marker size of the symbol for the target. Default is 75.
    """

```

2. Blind offset GMOS target with synthetic GAIA images

```python

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
    
    ***kwargs
    Remaining options are as similar to gmos_longslit_gaiasyn
    """
 ```

3. Two target GMOS acqusition with synthetic GAIA images

```python
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

    markersize: float, optional
    Markersize of User 1 and User 2 targets. Default is 75.
    
    markersizecentral: float, optional
    Markersize of the Target. Default is 25.
    
    """
```

Similar functions are available using the astroquery sky surveys. The available surveys are listed at https://astroquery.readthedocs.io/en/latest/skyview/skyview.html Default is set to DSS. The functions are

4. GMOS longslit target with astroquery skyview surveys

```python
    def gmos_longslit(name, ra1, dec1, pm_ra_cosdec=0.0, pm_dec=0.0, pa=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', skysurvey='DSS', frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75):


    """
    Function to create gmos finding charts from sky survey images, given the target coordinates, proper motions and epoch of observation.
    
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
    
    skysurvey: String, optional
    Skysurvey chosen from skyquery function. Options can be found at https://astroquery.readthedocs.io/en/latest/skyview/skyview.html. Defaults to DSS

    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.

    markersize: float, optional
    Markersize of the target. Default is 75.

    """
```

5. Blind offset GMOS target with astroquery skyview surveys

```python
    def gmos_blindoffset(name, ra1, dec1, rablind, decblind, pm_ra_cosdec=0.0, pm_dec=0.0, pa=0.0, slitwidth=1, pnum='Gemini GMOS Spectroscopy', skysurvey='DSS', frame='fk5', unit='deg', epoch='J2000', time_obs='J2020', corner='bottom left', markersize=75):


    """
    Function to create gmos finding charts from sky survey images, given the target coordinates, proper motions and epoch of observation.
    
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
    
    skysurvey: String, optional
    Skysurvey chosen from skyquery function. Options can be found at https://astroquery.readthedocs.io/en/latest/skyview/skyview.html. Defaults to DSS
    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.

    markersize: float, optional
    Markersize of the target. Default is 75.

    
    """
    
   ```
   
   6. GMOS two target acquisition using skyview
   
   ```python
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
    
    skysurvey: String, optional
    Skysurvey chosen from skyquery function. Options can be found at https://astroquery.readthedocs.io/en/latest/skyview/skyview.html. Defaults to DSS
    
    unit: string, or tuple, optional
    Units of the coordinates following astropy.Skycoord system. Default is decimal degrees.
    
    epoch: String, optional
    Epoch of the coordinates of the target. Default is set to J2000.
    
    time_obs: String, optional
    Epoch of when the obsevations will take place. Default is J2020.
    
    corner: String, optional
    Position of scalebar in the finding chart, following from aplpy convention
    Acceptable values are top right, top left, bottom right, bottom left, left, right, bottom or top.

    markersize: float, optional
    Markersize of User 1 and User 2 targets. Default is 75.
    
    markersizecentral: float, optional
    Markersize of the Target. Default is 25.

    
    """

```

Finally, if you dislike the finding charts, but would like to quickly generate GAIA synthetic images, you can use the function to create synthetic GAIA images. The field of view, and the FWHM cannot be changed. The FoV is 390 arcsec, while the FWHM is 0.3 arcsec. 

```python
    def make_gaia(coords, image=os.path.join('.','synthetic_gaia.fits'), epoch="J2020.5"):
    """
    Creates a synthetic image from GAIA DR2 photometry along the specified coordinates.

    Returns the synthetic image in fits format.

    Parameters:
    coords: Coordinates of centre of image. In astropy format. 
    
    image: File name of image. Default is synthetic_gaia.fits
    
    epoch: Epoch to translate image. Default is J2020.5
    
    """
```
