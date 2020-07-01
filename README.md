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

Creating a GMOS finding chart for a two target acqusition using synthetic GAIA DR2 images. You must enter the target name, and User 1 and User 2 coordinates. Optional arguements include  Program ID, slitwidth (default of 1 arcsec), and the position of the scalebar. The desired position angle, and the Base coordinates are calculated automatically, and saved to the log. Do not use the GAIA synthetic images if your target is a galaxy. 

    import gmosfc as gfc
    gfc.gmos_twotarget_gaiasyn("AL18_688+ALS18_689", '17h36m35.414s', '-33d30m12.7s', '17h36m36.11s', '-33d30m58.61s', slitwidth=2, pnum='GS-2020A-403', corner='bottom right')     
    ...
    
       
Creating a GMOS finding chart for a blind offset acqusition using DSS images. You must enter the target name, Base and User 1 coordinates. Additional arguements include the sky survey (deafults to DSS), which can be chosen from https://astroquery.readthedocs.io/en/latest/skyview/skyview.html  

    
    
    
For further information, please see documentation.     
    
Citing
------

If you have found this snippet useful, please give us a star.


Contact
----------

If you have issues with the code, please contact the author at vkalari@gemini.edu
