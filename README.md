gmosfc
==========

**Handy GMOS finding chart generator.**

gmosfc generates finding charts for Gemini-GMOS spectroscopy, using either archival sky survey images, or by generating synthetic images from GAIA DR2 photometry. 


Installing
----------

The following dependencies are required:

-  `Python <https://www.python.org/download/releases/3.0/>`__ 3.7.1 or later
-  `Numpy <http://www.numpy.org>`__ 1.11 or later
-  `Matplotlib <http://www.matplotlib.org>`__ 3.1.1 or later
-  `Astropy <http://www.astropy.org>`__ 4.0 or later
-  `Astroquery <https://astroquery.readthedocs.io/en/latest/>`__ 0.4 or later
-  `APLpy <https://aplpy.github.io/>`__ 2.0.2 or later
-  `photutils <https://photutils.readthedocs.io/en/stable/>`__ 0.7.2 or later

You can install APLpy and all its dependencies with:

    python3 -m pip install gmosfc
    


Examples
----------

Once installed, you can begin generating quick finding charts for a variety of GMOS modes.


*1. GMOS two target acqusition, using synthetic GAIA DR2 images*

Create a GMOS finding chart for a two target acqusition using synthetic GAIA DR2 images. You must enter the target name, and User 1 and User 2 coordinates in the astropy.SkyCoord format. The desired position angle, and the Base coordinates are calculated, and saved to the log.
Optional arguements include  Program ID, slitwidth (default of 1 arcsec), position of the scalebar, and SkyCoord options unit, frame, and epoch. Do not use the GAIA synthetic images if your target is a galaxy. 

    import gmosfc as gfc
    gfc.gmos_twotarget_gaiasyn("AL18_688+ALS18_689", '17h36m35.414s', '-33d30m12.7s', '17h36m36.11s', '-33d30m58.61s', slitwidth=2, pnum='GS-2020A-403', corner='bottom right')     
    ...
    Finding chart saved as AL18_688+ALS18_689_fc.jpg
    Cleaning up
    
The image is saved in the same directory, as the Target name fc.jpg    
![Image of ALS](https://github.com/astroquackers/gmosfc/blob/master/AL18_688%2BALS18_689_fc.jpg)

The utility of the GAIA synthetic images can be seen when you compare to the DSS finding chart created using the following command

    gfc.gmos_twotarget("AL18_688+ALS18_689", '17h36m35.414s', '-33d30m12.7s', '17h36m36.11s', '-33d30m58.61s', slitwidth=2, pnum='GS-2020A-403')
    
![Image of ALS](https://github.com/astroquackers/gmosfc/blob/master/AL18_688%2BALS18_689_dss.jpg)




*2. GMOS blind offset acqusition, using DSS images*

Creating a GMOS finding chart for a blind offset acqusition using DSS images. You must enter the target name, Base and User 1 coordinates. Additional arguements include the sky survey (deafults to DSS), which can be chosen from https://astroquery.readthedocs.io/en/latest/skyview/skyview.html  

    gfc.gmos_blindoffset('SDSSQuasar', '20h54m0.8016s', '-00d05m29.04s', '20h54m07.678s', '-00d05m31.92s')
    ...
    Finding chart saved as SDSSQuasar_fc.jpg    

We see now that the finding chart, and also that the default choice of scalebar position does not cover the compass.
![Image of QUASAR](https://github.com/astroquackers/gmosfc/blob/master/SDSSQuasar_fc.jpg)



*3. GMOS long slit target with position angle

A bright single target, with a desired position angle is entered as
    
    gfc.gmos_longslit('HD99', HD99', 001.492803, +44.739380, pa=123.4)
    
With the resulting finding chart
![HD99](https://github.com/astroquackers/gmosfc/blob/master/HD99_fc.jpg)
    

For further information, please see documentation.     
    
Citing
------

If you have found this snippet useful, please give us a star.


Contact
----------

If you have issues with the code, please contact the author at vkalari@gemini.edu
