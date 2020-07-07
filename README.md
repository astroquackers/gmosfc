gmosfc
==========

**Handy GMOS finding chart generator.**

gmosfc generates finding charts for Gemini-GMOS spectroscopy, using either archival sky survey images, or by generating synthetic images from GAIA DR2 photometry. gmosfc should be able to generate accurate finding charts for most GMOS spectroscopic modes or targets. If it is unable to do so, please contact us at vkalari@gemini.edu. 


Installing
----------

The following dependencies are required:

* [Python](https://www.python.org/download/releases/3.0/) 3.7.1 or later
* [Numpy](http://www.numpy.org) 1.11 or later
* [Matplotlib](http://www.matplotlib.org) 3.1.1 or later
* [Astropy](http://www.astropy.org) 4.0 or later
* [Astroquery](https://astroquery.readthedocs.io/en/latest/) 0.4 or later
* [APLpy](https://aplpy.github.io/) 2.0.2 or later
* [photutils](https://photutils.readthedocs.io/en/stable/) 0.7.2 or later

You can install APLpy and all its dependencies with:

    python3 -m pip install gmosfc
    


Examples
----------

Once installed, you can begin generating quick finding charts for a variety of GMOS modes.


*1. GMOS two target acqusition, using synthetic GAIA DR2 images*

Create a GMOS finding chart for a two target acqusition using synthetic GAIA DR2 images. You must enter the target name, and User 1 and User 2 coordinates in the astropy.SkyCoord format. The desired position angle, and the Base coordinates are calculated, and saved to the log.
Optional arguements include  Program ID, slitwidth (default of 1 arcsec), position of the scalebar, and SkyCoord options unit, frame, and epoch. 
*Do not use the GAIA synthetic images if your target is a galaxy.*

    import gmosfc as gfc
    gfc.gmos_twotarget_gaiasyn("Pismis24-17+Pismis24-16", 261.186378, -34.200760, 261.185225, -34.199680, slitwidth=2, pnum='GS-2020A-Q-403', markersize=25, markersizecentral=0)     
    ...
    Finding chart saved as Pismis24-17+Pismis24-16_fc.jpg
    Cleaning up
    
The image is saved in the same directory, as the Pismis24-17+Pismis24-16_fc.jpg    
![Image of ALS](https://github.com/astroquackers/gmosfc/blob/master/images/Pismis24-17%2BPismis24-16_fc_gaia.jpg)

The utility of the GAIA synthetic images can be seen when you compare to the DSS finding chart created using the following command (see the zoomed in comparison)

    gfc.gmos_twotarget("Pismis24-17+Pismis24-16", 261.186378, -34.200760, 261.185225, -34.199680, slitwidth=2, pnum='GS-2020A-Q-403', markersize=15, markersizecentral=0)
    
![Image of ALS](https://github.com/astroquackers/gmosfc/blob/master/images/Pismis24-17%2BPismis24-16_fc.jpg)

Looking at the zoomed in comparison of the GAIA synthetic and DSS finding charts, the image quality is apparent as the two targets are resolved (similar to what they would be in an 8m telescope image), as is the usefullness for the observer.  

![Image of Zoom](https://github.com/astroquackers/gmosfc/blob/master/images/comp.jpg)


*2. GMOS blind offset acqusition, using DSS images*

Creating a GMOS finding chart for a blind offset acqusition using DSS images. You must enter the target name, Base and User 1 coordinates. Additional arguements include the sky survey (deafults to DSS), which can be chosen from https://astroquery.readthedocs.io/en/latest/skyview/skyview.html  

    gfc.gmos_blindoffset('SDSSQuasar', '20h54m0.8016s', '-00d05m29.04s', '20h54m07.678s', '-00d05m31.92s')
    ...
    Finding chart saved as SDSSQuasar_fc.jpg    

We see now that the finding chart.
![Image of quasar](https://github.com/astroquackers/gmosfc/blob/master/images/SDSSQuasar_fc.jpg)


*3. GMOS long slit target with position angle*

A bright single target, with a desired position angle is entered as
    
    gfc.gmos_longslit('HD99', HD99', 001.492803, +44.739380, pa=123.4)
    
With the resulting finding chart
![HD99](https://github.com/astroquackers/gmosfc/blob/master/images/HD99_fc.jpg)
    

For further information, please see documentation.     

    
Citing
------

If you have found this package useful, please give us a star.
This code was funded partly by the Gemini-Conicyt grant 32RF180005, and the Gemini Observatory. 


Contact
----------

If you have issues with the code, please contact the author at vkalari@gemini.edu.


Disclaimer
----------

This is code is not officially endorsed software from the Gemini Observatory. The package is provided for convenience , and the author does not claim accuracy of the finding charts. It remains the responsibility of the PI to ensure that accurate finding charts are provided.



