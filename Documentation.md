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
    
We can know run functions avaiable in gmosfc as
 
    gfc.gmos_blindoffset()
    
The available functions for gfc using synthetic GAIA DR2 images are

